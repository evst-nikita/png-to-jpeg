#include "JPEGConverter.h"

// file-scope bit-buffer for writeBits/flushBitBuffer
static uint8_t g_bitBuffer = 0;  // здесь накапливаем очередные биты, сдвигая их в младшие позиции
static int   g_bitCount  = 0;    // сколько битов уже положили в g_bitBuffer


void JPEGConverter::loadImage(std::string filepath) {
    unsigned int w, h;
    unsigned int error = lodepng::decode(image, w, h, filepath);
    width = w;
    height = h;
    if (error) {
        std::cout << "Ошибка загрузки PNG: " << lodepng_error_text(error) << std::endl;
    }
}

void JPEGConverter::convert(string filepath) {
    loadImage(filepath);  // получение файла (png). Итоговый формат будет в виде вектора вида
    // [r, g, b, a, ...]

    blocksWidth = width / 16 + (width % 16 ? 1 : 0);  // считаем количество необходимых блоков
    blocksHeight = height / 16 + (height % 16 ? 1 : 0);

    Y3.resize(blocksWidth * blocksHeight * 16 * 16);  // это те же вектора но с большими числами
    Cb3.resize(blocksWidth * blocksHeight * 8 * 8);
    Cr3.resize(blocksWidth * blocksHeight * 8 * 8);
    {
        Y2.resize(blocksWidth * blocksHeight * 16 * 16);  // это те же вектора но с большими числами
        Cb2.resize(blocksWidth * blocksHeight * 8 * 8);
        Cr2.resize(blocksWidth * blocksHeight * 8 * 8);
        {
            Y.resize(blocksWidth * blocksHeight * 16 * 16);  // без сжатия
            Cb.resize(blocksWidth * blocksHeight * 8 * 8);  // сжатие вдвое
            Cr.resize(blocksWidth * blocksHeight * 8 * 8);

            makeBlocks();  // пакуем в блоки 8 на 8 и сразу приводим к диапазону -128 - 127

            DCTTransformBlocks();  // трансформируем блоки и применяем квантование
        }

        ZigZagTransform();  // трансформируем зиг загами и пакуем в векторы выше
    }

    unsigned int numMCU = blocksWidth * blocksHeight;
    unsigned int yBlocks = numMCU * 4;
    unsigned int cBlocks = numMCU;

    // Извлечение DC коэффициентов
    std::vector<int16_t> dcY(yBlocks), dcCb(cBlocks), dcCr(cBlocks);
    for (unsigned int i = 0; i < yBlocks; ++i) {
        dcY[i] = Y3[i * 64];
    }
    for (unsigned int i = 0; i < cBlocks; ++i) {
        dcCb[i] = Cb3[i * 64];
        dcCr[i] = Cr3[i * 64];
    }

    auto diffDCY  = encodeDiff(dcY);
    auto diffDCCb = encodeDiff(dcCb);
    auto diffDCCr = encodeDiff(dcCr);

    std::vector<uint8_t> dcSyms;
    dcSyms.reserve(diffDCY.size() + diffDCCb.size() + diffDCCr.size());
    for (auto d : diffDCY)  dcSyms.push_back(category(d));
    for (auto d : diffDCCb) dcSyms.push_back(category(d));
    for (auto d : diffDCCr) dcSyms.push_back(category(d));

    std::vector<uint8_t> acSyms;
    acSyms.reserve(yBlocks * 63 + 2 * cBlocks * 63);

    auto collectAC = [&](const std::vector<int16_t>& src, unsigned blocks){
        for (unsigned i = 0; i < blocks; ++i) {
            std::vector<int16_t> block(src.begin() + i*64 + 1, src.begin() + i*64 + 64);
            auto rle = countZeros(block);
            for (auto [run, val] : rle) {
                acSyms.push_back(acSymbolFromRunVal(run, val));
            }
        }
    };

    collectAC(Y3,  yBlocks);
    collectAC(Cb3, cBlocks);
    collectAC(Cr3, cBlocks);

    HuffmanEncoder dcEncoder, acEncoder;
    dcEncoder.build(dcSyms);
    acEncoder.build(acSyms);

    string newFilePath = filepath.substr(0, filepath.length() - 4) + "2.jpeg";
    ofstream out(newFilePath, std::ios::binary);
    if (out) {
        out.put((char) 0xFF);  // начальные байты файла
        out.put((char) 0xD8);
        writeApp0(out);
        writeQuantTables(out);  // запись таблиц квантования в файл
        writeSOF0(out);  // инфа по размерам, сжатию и таблицам

        {
            const auto& Li   = dcEncoder.getLi();
            const auto& Vals = dcEncoder.getVals();
            writeDHT(out, 0, 0, Li, Vals);
        }
        {
            const auto& Li   = acEncoder.getLi();
            const auto& Vals = acEncoder.getVals();
            writeDHT(out, 1, 0, Li, Vals);
        }
        writeSOS(out);
        BitWriter bw(out);

        int16_t prevY  = 0;
        int16_t prevCb = 0;
        int16_t prevCr = 0;

        if (yBlocks != 4 * cBlocks) {
            std::cerr << "Unexpected blocks layout: yBlocks=" << yBlocks
                      << " cBlocks=" << cBlocks << " (expected y=4*c)\n";
        }
        unsigned Nx = blocksWidth;              // число MCU по горизонтали
        unsigned widthBlocks = Nx * 2;          // число 8×8-блоков Y в строке
        unsigned numMCU = cBlocks;

        for (unsigned m = 0; m < numMCU; ++m) {
            unsigned yM = m / Nx, xM = m % Nx;
            unsigned base = (yM * 2) * widthBlocks + (xM * 2); // TL в сетке 8×8-блоков

            encodeBlock(bw, dcEncoder, acEncoder, Y3, base + 0,             prevY);
            encodeBlock(bw, dcEncoder, acEncoder, Y3, base + 1,             prevY);
            encodeBlock(bw, dcEncoder, acEncoder, Y3, base + widthBlocks + 0, prevY);
            encodeBlock(bw, dcEncoder, acEncoder, Y3, base + widthBlocks + 1, prevY);

            encodeBlock(bw, dcEncoder, acEncoder, Cb3, m, prevCb);
            encodeBlock(bw, dcEncoder, acEncoder, Cr3, m, prevCr);
        }

        bw.flush();

        out.put((char)0xFF); out.put((char)0xD9);
        out.close();
    } else {
        cout << "Error: can't save file";
    }
}

void JPEGConverter::encodeBlock(BitWriter& bw,
                               const HuffmanEncoder& dcEnc,
                               const HuffmanEncoder& acEnc,
                               const std::vector<int16_t>& src,
                               unsigned bi,
                               int16_t& prevDC)
{
    // --- DC ---
    int16_t DC   = src[bi*64 + 0];
    int16_t diff = DC - prevDC;           // предсказание от предыдущего DC той же компоненты
    prevDC       = DC;

    uint8_t sz   = category(diff);        // символ для DC — это size
    emitHuff(bw, dcEnc, sz);              // хаффман-код символа
    if (sz) {
        uint16_t amp = amplitudeBits(diff, sz); // доп.биты величины
        bw.putBits(amp, sz);
    }

    // --- AC ---
    // берём коэффициенты 1..63 этого блока
    std::vector<int16_t> ac(src.begin() + bi*64 + 1, src.begin() + bi*64 + 64);

    // ВАЖНО: твоя countZeros должна выдавать пары (run,val) и
    // ДОБАВИТЬ EOB (0,0), если в конце хвост нулей. Если не делает — см. заметку ниже.
    auto rle = JPEGConverter::countZeros(ac);

    for (const auto& p : rle) {
        int8_t  run = p.first;
        int16_t val = p.second;

        if (val == 0) {
            // либо ZRL (15,0) => 0xF0, либо EOB => 0x00
            emitHuff(bw, acEnc, (run == 15) ? 0xF0 : 0x00);
        } else {
            uint8_t s  = category(val);
            uint8_t sym = (uint8_t)((run << 4) | (s & 0x0F));
            emitHuff(bw, acEnc, sym);
            bw.putBits(amplitudeBits(val, s), s);
        }
    }
}


void JPEGConverter::writeSOS(std::ofstream& out) {
    out.put((char)0xFF); out.put((char)0xDA);
    const uint8_t Ns = 3;
    const uint16_t L = 2 + 1 + Ns*2 + 3; // 2 + Ns + 2*Ns + 3
    out.put((char)(L >> 8)); out.put((char)L);
    out.put((char)Ns);

    // id=1(Y),2(Cb),3(Cr); Td=0, Ta=0
    struct { uint8_t id, td, ta; } cs[3] = {{1,0,0},{2,0,0},{3,0,0}};
    for (auto c : cs) {
        out.put((char)c.id);
        out.put((char)((c.td << 4) | (c.ta & 0x0F)));
    }

    out.put((char)0x00); out.put((char)0x3F); out.put((char)0x00);
}


void JPEGConverter::writeDHT(std::ofstream& out, uint8_t tc, uint8_t th,
                             const std::vector<uint16_t>& Li,     // Li[1..16]
                             const std::vector<uint8_t>&  Vals)   // huffval[]
{
    auto sumLi = size_t{0};
    for (int i = 1; i <= 16; ++i) {
        if (i < (int)Li.size()) sumLi += Li[i];
    }
    if (sumLi != Vals.size()) {
        std::cerr << "[DHT ERROR] sum(Li)=" << sumLi
                  << " != Vals.size()=" << Vals.size() << '\n';
    }
    size_t nwrite = Vals.size();


    // Ls = 2 (поле длины) + 1 (Tc/Th) + 16 (Li) + nwrite
    uint16_t Ls = uint16_t(2 + 1 + 16 + nwrite);

    out.put((char)0xFF); out.put((char)0xC4);
    out.put((char)(Ls >> 8)); out.put((char)Ls);
    out.put((char)((tc << 4) | (th & 0x0F)));

    for (int i = 1; i <= 16; ++i) {
        uint16_t cnt = (i < (int)Li.size()) ? Li[i] : 0;
        out.put((char)std::min<uint16_t>(cnt, 255));
    }

    // Пишем ровно nwrite значений, не больше
    out.write(reinterpret_cast<const char*>(Vals.data()), std::streamsize(nwrite));

    // Диагностика на всякий случай
    if (sumLi != Vals.size()) {
        std::cerr << "[DHT warn] sum(Li)=" << sumLi
                  << " != Vals.size()=" << Vals.size()
                  << " -> wrote " << nwrite << " bytes\n";
    }
}


inline void JPEGConverter::emitHuff(BitWriter& bw, const HuffmanEncoder& enc, uint8_t sym) {
    uint32_t code; int len;
    if (!enc.getCode(sym, code, len)) {
        return;
    }
    bw.putBits(code, len); // MSB-first
}

void JPEGConverter::encodeDC(BitWriter& bw, const HuffmanEncoder& dcEnc, int16_t diff) {
    uint8_t sz = category(diff);
    emitHuff(bw, dcEnc, sz);
    if (sz) {
        uint16_t bits = amplitudeBits(diff, sz);
        bw.putBits(bits, sz);
    }
}

void JPEGConverter::encodeAC(BitWriter& bw, const HuffmanEncoder& acEnc,
                            const std::vector<int16_t>& src, unsigned bi,
                            std::vector<std::pair<int8_t,int16_t>> (*countZeros)(const std::vector<int16_t>&))
{
    std::vector<int16_t> ac(src.begin() + bi*64 + 1, src.begin() + bi*64 + 64);
    auto rle = countZeros(ac);

    for (const auto& p : rle) {
        int8_t  run = p.first;
        int16_t val = p.second;

        if (val == 0) {
            uint8_t sym = (run == 15) ? 0xF0 : 0x00;
            emitHuff(bw, acEnc, sym);
        } else {
            uint8_t sz  = category(val);
            uint8_t sym = (uint8_t)((run << 4) | (sz & 0x0F));
            emitHuff(bw, acEnc, sym);
            uint16_t bits = amplitudeBits(val, sz);
            bw.putBits(bits, sz);
        }
    }
}

uint16_t JPEGConverter::amplitudeBits(int16_t v, uint8_t sz) {
    if (sz == 0) return 0;
    if (v >= 0) return (uint16_t)v;
    uint16_t mask = (1u << sz) - 1u;
    return (~(uint16_t)(-v)) & mask; // ones' complement
}

uint8_t JPEGConverter::category(int16_t v) {
    if (v == 0) return 0;
    uint16_t a = (v < 0) ? -v : v;
    uint8_t n = 0;
    while (a) { a >>= 1; ++n; }
    return n;
}

uint8_t JPEGConverter::acSymbolFromRunVal(int8_t run, int16_t val) {
    if (val == 0) {
        return (run == 15) ? 0xF0 : 0x00;
    } else {
        uint8_t sz = category(val);
        return (uint8_t)((run << 4) | (sz & 0x0F));
    }
}

void JPEGConverter::makeBlocks() {
    // Блок работы с яркостью
    unsigned int YIndex = 0;
    for (unsigned int i = 0; i < height; ++i) {
        for (unsigned int j = 0; j < width; ++j) {
            unsigned int index = 4 * (i * width + j);  // тут 4 т.к. 4 канала
            uint8_t r = image[index];
            uint8_t g = image[index + 1];
            uint8_t b = image[index + 2];
            int color = (int) round(0.299 * r + 0.587 * g + 0.114 * b);
            if (color > 255) color = 255;
            Y[YIndex++] = (int8_t) (color - 128);
        }
        for (unsigned int j = width; j < blocksWidth * 16; ++j) {  // добивание остальных клеток по ширине
            Y[YIndex++] = 0;
        }
    }
    for (unsigned int i = height; i < blocksHeight * 16; ++i) {  // добивание нижних рядов
        for (unsigned int j = 0; j < blocksWidth * 16; ++j) {
            Y[YIndex++] = 0;
        }
    }

    vector<int8_t> CbExtra, CrExtra;
    CbExtra.resize(blocksWidth * blocksHeight * 16 * 16);
    CrExtra.resize(blocksWidth * blocksHeight * 16 * 16);
    unsigned int CIndex = 0;
    for (unsigned int i = 0; i < height; ++i) {
        for (unsigned int j = 0; j < width; ++j) {
            unsigned int index = 4 * (i * width + j);  // тут 4 т.к. 4 канала
            uint8_t r = image[index];
            uint8_t g = image[index + 1];
            uint8_t b = image[index + 2];
            int CbColor = (int) round(128 - 0.168736 * r - 0.331264 * g + 0.5 * b);
            if (CbColor > 255) CbColor = 255;
            if (CbColor < 0) CbColor = 0;
            int CrColor = (int) round(128 + 0.5 * r - 0.418688 * g - 0.081312 * b);
            if (CrColor > 255) CrColor = 255;
            if (CrColor < 0) CrColor = 0;
            CbExtra[CIndex] = (int8_t) (CbColor - 128);
            CrExtra[CIndex++] = (int8_t) (CrColor - 128);
        }
        for (unsigned int j = width; j < blocksWidth * 16; ++j) {  // дозаполнение остальных клеток
            CbExtra[CIndex] = 0;
            CrExtra[CIndex++] = 0;
        }
    }
    for (unsigned int i = height; i < blocksHeight * 16; ++i) {
        for (unsigned int j = 0; j < blocksWidth * 16; ++j) {
            CbExtra[CIndex] = 0;
            CrExtra[CIndex++] = 0;
        }
    }

    // Теперь сжатие Cb и Cr
    CIndex = 0;
    for (unsigned int i = 0; i < blocksHeight * 16; i += 2) {
        for (unsigned int j = 0; j < blocksWidth * 16; j += 2) {
            unsigned int index = i * blocksWidth * 16 + j;
            unsigned int index2 = (i + 1) * blocksWidth * 16 + j;
            Cb[CIndex] = (int8_t) round((double) (CbExtra[index] + CbExtra[index + 1] + CbExtra[index2] + CbExtra[index2 + 1]) / 4);
            Cr[CIndex++] = (int8_t) round((double) (CrExtra[index] + CrExtra[index + 1] + CrExtra[index2] + CrExtra[index2 + 1]) / 4);
        }
    }
}

void JPEGConverter::printY2() {
    for (int i = 0; i < blocksHeight * 8; ++i) {
        for (int j = 0; j < blocksWidth * 8; ++j) {
            cout << (int) Cr2[i * blocksWidth * 8 + j] << ' ';
        }
        cout << endl;
    }
}

void JPEGConverter::DCTTransformBlocks() {
    vector<int8_t> buffer(64);
    vector<int8_t> buffer2(64);

    for (unsigned int i = 0; i < blocksHeight * 16; i += 8) {
        for (unsigned int j = 0; j < blocksWidth * 16; j += 8) {
            for (int y = 0; y < 8; ++y) {  // пакуем блок в массив из 64 элементов
                for (int x = 0; x < 8; ++x) {
                    buffer[y * 8 + x] = Y[(i + y) * blocksWidth * 16 + j + x];
                }
            }
            vector<int16_t> newBuffer = DCTTransform(buffer, true);
            for (int y = 0; y < 8; ++y) {  // распаковываем блок в массив из 64 элементов
                for (int x = 0; x < 8; ++x) {
                    Y2[(i + y) * blocksWidth * 16 + j + x] = newBuffer[y * 8 + x];
                }
            }
        }
    }
    for (unsigned int i = 0; i < blocksHeight * 8; i += 8) {
        for (unsigned int j = 0; j < blocksWidth * 8; j += 8) {
            for (int y = 0; y < 8; ++y) {  // пакуем блок в массив из 64 элементов
                for (int x = 0; x < 8; ++x) {
                    buffer[y * 8 + x] = Cb[(i + y) * blocksWidth * 8 + j + x];
                    buffer2[y * 8 + x] = Cr[(i + y) * blocksWidth * 8 + j + x];
                }
            }
            vector<int16_t> newBuffer = DCTTransform(buffer, false);
            vector<int16_t> newBuffer2 = DCTTransform(buffer2, false);
            for (int y = 0; y < 8; ++y) {  // распаковываем блок в массив из 64 элементов
                for (int x = 0; x < 8; ++x) {
                    Cb2[(i + y) * blocksWidth * 8 + j + x] = newBuffer[y * 8 + x];
                    Cr2[(i + y) * blocksWidth * 8 + j + x] = newBuffer2[y * 8 + x];
                }
            }
        }
    }
}

double JPEGConverter::alpha(int omega) {
    if (omega == 0)
        return 1.0 / sqrt(2);
    else
        return 1;
}

vector<int16_t> JPEGConverter::DCTTransform(vector<int8_t> &buffer, bool isLum) {
    vector<double> DCTbuffer(64);
    vector<int16_t> result(64);

    for (int v = 0; v < 8; ++v) {  // v - вертикаль
        for (int u = 0; u < 8; ++u) {  // u - горизонталь
            double sum = 0;
            for (int y = 0; y < 8; ++y) {
                for (int x = 0; x < 8; ++x) {
                    sum += (double) buffer[y * 8 + x] * cosTable[x * 8 + u] * cosTable[y * 8 + v];
                }
            }
            DCTbuffer[v * 8 + u] = sum * 0.25 * alpha(u) * alpha(v);
        }
    }

    if (isLum) {  // фактически выбор таблицы квантования
        for (int y = 0; y < 8; ++y) {
            for (int x = 0; x < 8; ++x) {
                result[y * 8 + x] = (int16_t) round(DCTbuffer[y * 8 + x] / quantLuminance[y][x]);
            }
        }
    } else {
        for (int y = 0; y < 8; ++y) {
            for (int x = 0; x < 8; ++x) {
                result[y * 8 + x] = (int16_t) round(DCTbuffer[y * 8 + x] / quantChrominance[y][x]);
            }
        }
    }
    return result;
}

JPEGConverter::JPEGConverter() {
    constexpr double PI = 3.14159265358979323846;
    cosTable.resize(64);
    for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
            cosTable[a * 8 + b] = cos((double) (2 * a + 1) * b * PI / 16);
        }
    }
}

void JPEGConverter::ZigZagTransform() {
    vector<int16_t> buffer(64);
    vector<int16_t> buffer2(64);
    int blockIdx = 0;
    for (unsigned int i = 0; i < blocksHeight * 16; i += 8) {
        for (unsigned int j = 0; j < blocksWidth * 16; j += 8) {
            for (int y = 0; y < 8; ++y) {  // пакуем блок в массив из 64 элементов
                for (int x = 0; x < 8; ++x) {
                    buffer[y * 8 + x] = Y2[(i + y) * blocksWidth * 16 + j + x];
                }
            }
            for (int k = 0; k < 64; ++k) {
                Y3[blockIdx * 64 + k] = buffer[zigZagOrder[k]];
            }
            blockIdx += 1;
        }
    }

    blockIdx = 0;

    for (unsigned int i = 0; i < blocksHeight * 8; i += 8) {
        for (unsigned int j = 0; j < blocksWidth * 8; j += 8) {
            for (int y = 0; y < 8; ++y) {  // пакуем блок в массив из 64 элементов
                for (int x = 0; x < 8; ++x) {
                    buffer[y * 8 + x] = Cb2[(i + y) * blocksWidth * 8 + j + x];
                    buffer2[y * 8 + x] = Cr2[(i + y) * blocksWidth * 8 + j + x];
                }
            }
            // распаковываем блок в массив из 64 элементов
            for (int k = 0; k < 64; ++k) {
                Cb3[blockIdx * 64 + k] = buffer[zigZagOrder[k]];
                Cr3[blockIdx * 64 + k] = buffer2[zigZagOrder[k]];
            }
            blockIdx += 1;
        }
    }
}

vector<pair<int8_t, int16_t>> JPEGConverter::countZeros(const vector<int16_t>& input) {
    vector<pair<int8_t, int16_t>> res;
    int8_t zeros = 0;
    for (int16_t num : input) {
        if (num != 0 || zeros == 15) {
            res.emplace_back(zeros, num);
            zeros = 0;
        } else {
            ++zeros;
        }
    }
    if (zeros != 0)
        res.emplace_back(0, 0);
    return res;
}

vector<int16_t> JPEGConverter::encodeDiff(const std::vector<int16_t>& dc) {
    std::vector<int16_t> diff;
    diff.reserve(dc.size());
    int16_t prev = 0;
    for (auto v : dc) {
        diff.push_back(v - prev);
        prev = v;
    }
    return diff;
}

void JPEGConverter::writeQuantTables(ofstream &out) {
    out.put((char) 0xFF);
    out.put((char) 0xDB);
    // общая длина будет 2 + 2 * (64 + 1) = 132 байта
    out.put((char) 0x00);
    out.put((char) 0x84);
    // Далее вставляем первую таблицу
    out.put((char) 0x00);
    for (int i = 0; i < 64; ++i)
        out.put((char) quantLuminance[zigZagOrder[i] / 8][zigZagOrder[i] % 8]);
    // вставляем вторую таблицу
    out.put((char) 0x01);
    for (int i = 0; i < 64; ++i)
        out.put((char) quantChrominance[zigZagOrder[i] / 8][zigZagOrder[i] % 8]);
}

void JPEGConverter::writeApp0(ofstream &out) {
    uint8_t app0[] = {
            0xFF, 0xE0,       // маркер APP0
            0x00, 0x10,       // длина = 16
            0x4A,0x46,0x49,0x46,0x00, // "JFIF\0"
            0x01,0x00,        // версия 1.00
            0x00,             // Units = 0
            0x00,0x01,        // Xdensity = 1
            0x00,0x01,        // Ydensity = 1
            0x00,             // Xthumbnail = 0
            0x00              // Ythumbnail = 0
    };
    for (auto b : app0) out.put((char) b);
}

void JPEGConverter::writeSOF0(ofstream &out) {
    // Маркер SOF0
    out.put(char(0xFF)); out.put(char(0xC0));

    // Вычисляем длину: 2 + 1 + 2 + 2 + 1 + numComps*3
    constexpr int numComps = 3;
    uint16_t len = 2 + 1 + 2 + 2 + 1 + numComps * 3;
    out.put((char) (len >> 8)); out.put((char) len);

    // Precision
    out.put(char(0x08));

    // Высота и ширина (big-endian)
    out.put((char) (height >> 8)); out.put((char) height);
    out.put((char) (width >> 8)); out.put((char) width);

    // Число компонентов
    out.put((char) numComps);

    // Для каждого компонента: ID, sampling-факторы, Tq
    struct Comp { uint8_t id, h, v, tq; };
    Comp comps[numComps] = {
            {1, 2, 2, 0}, // Y, Tq=0
            {2, 1, 1, 1}, // Cb, Tq=1
            {3, 1, 1, 1}  // Cr, Tq=1
    };
    for (const auto &c : comps) {
        out.put((char) c.id);
        out.put((char) ((c.h << 4) | c.v));
        out.put((char) c.tq);
    }
}

BitWriter::BitWriter(std::ostream& out) : out_(out) {}

void BitWriter::putByte(uint8_t b) {
        out_.put((char)b);
        if (b == 0xFF) out_.put((char)0x00);
    }

void BitWriter::putBits(uint32_t bits, int n) {
    for (int i = n-1; i >= 0; --i) {
        buf_ = (buf_ << 1) | ((bits >> i) & 1u);
        if (++cnt_ == 8) { putByte(buf_); buf_ = 0; cnt_ = 0; }
    }
}

void BitWriter::flush() {
    if (cnt_ > 0) { buf_ <<= (8 - cnt_); putByte(buf_); buf_ = 0; cnt_ = 0; }
}
