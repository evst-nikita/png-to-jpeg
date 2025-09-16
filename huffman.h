#pragma once
#include <vector>
#include <queue>
#include <unordered_map>
#include <cstdint>
#include <algorithm>

// Класс строит канонические хаффман-коды по распределению байтовых символов (0..255).
// ВАЖНО для JPEG: длины кодов ограничиваются 16 битами по процедуре Annex K (Figure K.3).
class HuffmanEncoder {
public:
    using Symbol = uint8_t;

    // Построить таблицу по массиву символов (байт 0..255)
    void build(const std::vector<Symbol>& data) {
        clear_();

        // 1) Частоты
        std::vector<size_t> freq(256, 0);
        for (auto s : data) freq[s]++;

        // Пустой поток — пустые таблицы
        bool any = false;
        for (auto f : freq) { if (f) { any = true; break; } }
        if (!any) {
            Li_.assign(17, 0); Vals_.clear(); codes_.clear(); codeLengths_.clear();
            return;
        }

        // 2) Построить дерево Хаффмана (минимальная куча по (freq, symbol))
        struct Node {
            size_t f; int s; Node* l; Node* r;
        };
        auto cmp = [](Node* a, Node* b){
            if (a->f != b->f) return a->f > b->f;
            return a->s > b->s;
        };
        std::priority_queue<Node*, std::vector<Node*>, decltype(cmp)> pq(cmp);
        std::vector<Node*> nodes; nodes.reserve(512);

        for (int s = 0; s < 256; ++s) if (freq[s]) {
                Node* n = new Node{freq[s], s, nullptr, nullptr};
                pq.push(n); nodes.push_back(n);
            }
        // Спец-кейс: один символ -> фиктивный родитель, чтобы длина стала 1
        if (pq.size() == 1) {
            Node* only = pq.top(); pq.pop();
            Node* parent = new Node{only->f, only->s, only, nullptr};
            pq.push(parent); nodes.push_back(parent);
        }
        while (pq.size() > 1) {
            Node* a = pq.top(); pq.pop();
            Node* b = pq.top(); pq.pop();
            int sym = (a->s < b->s) ? a->s : b->s; // стабилизация по символу
            Node* p = new Node{a->f + b->f, sym, a, b};
            pq.push(p); nodes.push_back(p);
        }
        Node* root = pq.top();

        // 3) Исходные длины кодов (глубина листа)
        assignLengths_(root, 0);
        // Для single-symbol глубина 0 -> сделать 1
        for (auto& kv : codeLengths_) if (kv.second == 0) kv.second = 1;

        // 4) Подсчёт количества кодов каждой длины
        //   bits[L] = сколько символов с длиной L. Длины могут быть >16 -> отложим ограничение.
        std::vector<int> bits(33, 0); // индексы 1..32 (запас)
        int maxLen = 0;
        for (auto& kv : codeLengths_) {
            int L = kv.second;
            if (L < 1) L = 1;
            if (L > 32) L = 32; // страховка, до перераспределения
            bits[L]++;
            if (L > maxLen) maxLen = L;
        }

        // 5) Ограничить длины до 16 по Annex K (Figure K.3)
        // Идея: пока есть коды длиной >16, «подрезаем» их, продвигая часть кодов вверх,
        // сохраняя префиксные свойства (см. libjpeg: jpeg_gen_optimal_table).
        for (int i = maxLen; i > 16; --i) {
            while (bits[i] > 0) {
                int j = i - 2;
                while (j > 0 && bits[j] == 0) --j;

                // ЗАЩИТА от underflow/некорректных переносов
                if (bits[i] < 2 || j <= 0 || bits[j] == 0) break;

                // убрать два самых длинных кода
                bits[i]      -= 2;
                // один код переносим на длину i-1
                bits[i - 1]  += 1;
                // из длины j «повышаем» префикс
                bits[j]      -= 1;
                bits[j + 1]  += 2;
            }
        }


        // 6) Сформировать Li[1..16] для DHT
        Li_.assign(17, 0);
        for (int L = 1; L <= 16; ++L) Li_[L] = bits[L];

        // 6a) Избежать «all-ones» кодов (libjpeg: JERR_BAD_HUFF_TABLE)
        auto avoid_all_ones = [&](){
            std::vector<int> li(17, 0);
            for (int L = 1; L <= 16; ++L) li[L] = (int)Li_[L];

            auto has_bad = [&](int& badL)->bool{
                unsigned code = 0;
                for (int L = 1; L <= 16; ++L) {
                    code = (code + li[L-1]) << 1;             // стартовый код длины L
                    if (li[L] > 0) {
                        unsigned end = code + (unsigned)li[L] - 1u; // последний код длины L
                        if (end == ((1u << L) - 1u)) {        // все биты 1 → плохо
                            badL = L;
                            return true;
                        }
                    }
                }
                return false;
            };

            int badL = 0, guard = 0;
            while (has_bad(badL) && guard++ < 64) {
                // найдём ближайший более короткий уровень с ненулевым количеством
                int k = badL - 1;
                while (k >= 1 && li[k] == 0) --k;
                if (k < 1) break;              // не нашли — выходим

                // переносим ОДИН символ с уровня k на самый глубокий уровень (16),
                // это создаёт "запас" и уводит последний код от границы all-ones
                li[k] -= 1;
                li[16] += 1;
            }

            for (int L = 1; L <= 16; ++L) Li_[L] = (uint16_t)li[L];
        };
        avoid_all_ones();


        // 7) Распределить символы по новым длинам:
        //    сортируем символы по (исходная длина, затем по значению символа),
        //    затем по очереди выдаём их в группы длины 1..16 согласно Li.
        std::vector<Symbol> symbols;
        symbols.reserve(codeLengths_.size());
        for (auto& kv : codeLengths_) symbols.push_back((Symbol)kv.first);
        std::sort(symbols.begin(), symbols.end(), [&](Symbol a, Symbol b){
            int la = codeLengths_[a], lb = codeLengths_[b];
            if (la != lb) return la < lb;
            return a < b;
        });

        std::vector<std::vector<Symbol>> byLen(17);
        {
            size_t pos = 0;
            for (int L = 1; L <= 16; ++L) {
                int take = Li_[L];
                while (take > 0 && pos < symbols.size()) {
                    byLen[L].push_back(symbols[pos]);
                    codeLengths_[symbols[pos]] = L; // обновили длину на «уложенную»
                    ++pos; --take;
                }
            }
            // На всякий случай: отсортируем внутри каждой длины по возрастанию символа (как в DHT)
            for (int L = 1; L <= 16; ++L) {
                std::sort(byLen[L].begin(), byLen[L].end());
            }
            // Синхронизируем Li_ с фактическими размерами групп
            for (int L = 1; L <= 16; ++L) {
                Li_[L] = static_cast<uint16_t>(byLen[L].size());
            }

        }

        // 8) Сформировать Vals_ в порядке длин 1..16
        Vals_.clear(); Vals_.reserve(symbols.size());
        for (int L = 1; L <= 16; ++L) {
            Vals_.insert(Vals_.end(), byLen[L].begin(), byLen[L].end());
        }

        // 9) Назначить канонические коды
        // Начальные коды по длинам (см. Annex C / libjpeg)
        codes_.clear();
        std::vector<uint16_t> nextCode(17, 0);
        uint16_t code = 0;
        // Li_[0] считаем 0
        for (int L = 1; L <= 16; ++L) {
            code = (code + Li_[L - 1]) << 1;
            nextCode[L] = code;
        }
        // Идём по Vals_ в порядке длин (это именно порядок DHT)
        size_t idx = 0;
        for (int L = 1; L <= 16; ++L) {
            for (int n = 0; n < Li_[L]; ++n, ++idx) {
                Symbol s = Vals_[idx];
                codes_[s] = { nextCode[L], L };
                nextCode[L]++;
            }
        }

        // 10) Очистка узлов дерева
        for (auto* n : nodes) delete n;
    }

    // Получить (код, длину) для символа. Возвращает false, если символа нет в таблице.
    bool getCode(Symbol s, uint32_t& code, int& len) const {
        auto it = codes_.find(s);
        if (it == codes_.end()) return false;
        code = it->second.first; len = it->second.second;
        return true;
    }

    // Удобная «склейка» битов в вектор<bool> (для тестов)
    std::vector<bool> encode(const std::vector<Symbol>& data) const {
        std::vector<bool> out;
        for (auto s : data) {
            uint32_t c; int l;
            if (getCode(s, c, l)) {
                for (int i = l - 1; i >= 0; --i) out.push_back((c >> i) & 1u);
            }
        }
        return out;
    }

    // Для записи DHT
    const std::vector<uint16_t>& getLi()   const { return Li_; }   // Li_[1..16]
    const std::vector<Symbol>&   getVals() const { return Vals_; } // символы в порядке длин

private:
    // длины исходного дерева (ещё до укладки ≤16)
    std::unordered_map<Symbol, int> codeLengths_;               // длины кодов
    std::unordered_map<Symbol, std::pair<uint32_t,int>> codes_; // symbol -> (code,len)
    std::vector<uint16_t> Li_;   // 17 элементов, индексы 1..16
    std::vector<Symbol>   Vals_; // последовательность символов по длинам

    void clear_() {
        codeLengths_.clear(); codes_.clear(); Li_.clear(); Vals_.clear();
    }

    // Рекурсивно назначаем длины кодов (глубина листа)
    void assignLengths_(void* vnode, int d) {
        struct Node { size_t f; int s; Node* l; Node* r; };
        Node* node = static_cast<Node*>(vnode);
        if (!node) return;
        if (!node->l && !node->r) {
            codeLengths_[(Symbol)node->s] = d;
            return;
        }
        assignLengths_(node->l, d + 1);
        assignLengths_(node->r, d + 1);
    }

    // NEW: инициализация из готовых Li/Vals (Annex K)
    void initFromLiVals(const std::vector<uint16_t>& Li, const std::vector<uint8_t>& Vals) {
        // сохраняем
        Li_.assign(17, 0);
        for (int i = 1; i <= 16 && i < (int)Li.size(); ++i) Li_[i] = Li[i];
        Vals_ = Vals;

        // построить канонические коды
        codes_.clear();
        std::vector<uint16_t> nextCode(17, 0);
        uint16_t code = 0;
        for (int L = 1; L <= 16; ++L) {
            code = (code + Li_[L-1]) << 1;   // Li_[0] = 0
            nextCode[L] = code;
        }
        size_t pos = 0;
        for (int L = 1; L <= 16; ++L) {
            for (int n = 0; n < Li_[L]; ++n, ++pos) {
                Symbol s = Vals_[pos];
                codes_[s] = { nextCode[L], L };
                nextCode[L]++;
            }
        }
    }

};
