#pragma once
#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include "lodepng.h"
#include <string>
#include <utility>
#include "huffman.h"
#include <fstream>
#include <numeric>

using std::cos;
using std::cout;
using std::endl;
using std::ofstream;
using std::ostream;
using std::pair;
using std::string;
using std::sqrt;
using std::vector;

class BitWriter;

class JPEGConverter {
private:
    vector<uint8_t> image;
    uint16_t width, height;
    unsigned int blocksWidth, blocksHeight;
    vector<int8_t> Y, Cb, Cr;
    vector<int16_t> Y2, Cb2, Cr2;
    vector<int16_t> Y3, Cb3, Cr3;
    vector<double> cosTable;
    const uint8_t quantLuminance[8][8] = {
            { 16, 11, 10, 16,  24,  40,  51,  61 },
            { 12, 12, 14, 19,  26,  58,  60,  55 },
            { 14, 13, 16, 24,  40,  57,  69,  56 },
            { 14, 17, 22, 29,  51,  87,  80,  62 },
            { 18, 22, 37, 56,  68, 109, 103,  77 },
            { 24, 35, 55, 64,  81, 104, 113,  92 },
            { 49, 64, 78, 87, 103, 121, 120, 101 },
            { 72, 92, 95, 98, 112, 100, 103,  99 }
    };
    const uint8_t quantChrominance[8][8] = {
            { 17, 18, 24, 47, 99,  99,  99,  99 },
            { 18, 21, 26, 66, 99,  99,  99,  99 },
            { 24, 26, 56, 99, 99,  99,  99,  99 },
            { 47, 66, 99, 99, 99,  99,  99,  99 },
            { 99, 99, 99, 99, 99,  99,  99,  99 },
            { 99, 99, 99, 99, 99,  99,  99,  99 },
            { 99, 99, 99, 99, 99,  99,  99,  99 },
            { 99, 99, 99, 99, 99,  99,  99,  99 }
    };
    const uint8_t zigZagOrder[64] = {
            0,  1,  8, 16,  9,  2,  3, 10,
            17, 24, 32, 25, 18, 11,  4,  5,
            12, 19, 26, 33, 40, 48, 41, 34,
            27, 20, 13,  6,  7, 14, 21, 28,
            35, 42, 49, 56, 57, 50, 43, 36,
            29, 22, 15, 23, 30, 37, 44, 51,
            58, 59, 52, 45, 38, 31, 39, 46,
            53, 60, 61, 54, 47, 55, 62, 63
    };

public:
    JPEGConverter();

    void loadImage(string filepath);
    void convert(string filepath);
    void makeBlocks();
    void DCTTransformBlocks();
    void ZigZagTransform();
    void writeQuantTables(ofstream &out);
    void writeApp0(ofstream &out);
    void writeSOF0(ofstream &out);
    void writeDHT(ofstream &out, uint8_t tc, uint8_t th, const vector<uint16_t> &Li, const vector<uint8_t> &Vals);
    void writeSOS(std::ofstream& out);
    void encodeBlock(BitWriter &bw, const HuffmanEncoder &dcEnc, const HuffmanEncoder &acEnc, const vector<int16_t> &src,
                     unsigned int bi, int16_t &prevDC);

    vector<int16_t> DCTTransform(vector<int8_t>& buffer, bool isLum);
    static double alpha(int omega);
    uint8_t category(int16_t v);
    uint8_t acSymbolFromRunVal(int8_t run, int16_t val);
    void emitHuff(BitWriter& bw, const HuffmanEncoder& enc, uint8_t sym);
    void encodeDC(BitWriter& bw, const HuffmanEncoder& dcEnc, int16_t diff);
    void encodeAC(BitWriter& bw, const HuffmanEncoder& acEnc,
                                 const std::vector<int16_t>& src, unsigned bi,
                                 std::vector<std::pair<int8_t,int16_t>> (*countZeros)(const std::vector<int16_t>&));
    static vector<pair<int8_t, int16_t>> countZeros(const vector<int16_t>& input);
    vector<int16_t> encodeDiff(const std::vector<int16_t>& dc);
    uint16_t amplitudeBits(int16_t v, uint8_t sz);

    void printY2();
};

class BitWriter {
public:
    BitWriter(ostream& out);
    void putByte(uint8_t b);
    void putBits(uint32_t bits, int n);
    void flush();
private:
    ostream& out_;
    uint8_t buf_ = 0;
    int     cnt_ = 0;
};

