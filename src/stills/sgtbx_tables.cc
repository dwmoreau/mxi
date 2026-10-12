// cctbx's space-group reference tables, copied verbatim from cctbx_project
// cctbx/sgtbx/reference_settings/{hall_symbol_table,matrix_group_code_table,
// normalizer}.cpp and lattice_symmetry.h by
// claude/scripts/rewrite/r8a_sgtbx_tables.py in stills_process_profiling. Do
// not edit by hand.
//
// cctbx Copyright (c) 2006 - 2026, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to
// receipt of any required approvals from the U.S. Dept. of Energy). All
// rights reserved. BSD-3-Clause-LBNL; see third_party/cctbx/LICENSE.txt.

#include "stills/sgtbx.hh"

namespace mxi::stills::sgtbx {

const char *hall_symbol_table(int number) {
  static const char *const table[] = {
      nullptr,
      " P 1",            // 001
      "-P 1",            // 002
      " P 2y",           // 003
      " P 2yb",          // 004
      " C 2y",           // 005
      " P -2y",          // 006
      " P -2yc",         // 007
      " C -2y",          // 008
      " C -2yc",         // 009
      "-P 2y",           // 010
      "-P 2yb",          // 011
      "-C 2y",           // 012
      "-P 2yc",          // 013
      "-P 2ybc",         // 014
      "-C 2yc",          // 015
      " P 2 2",          // 016
      " P 2c 2",         // 017
      " P 2 2ab",        // 018
      " P 2ac 2ab",      // 019
      " C 2c 2",         // 020
      " C 2 2",          // 021
      " F 2 2",          // 022
      " I 2 2",          // 023
      " I 2b 2c",        // 024
      " P 2 -2",         // 025
      " P 2c -2",        // 026
      " P 2 -2c",        // 027
      " P 2 -2a",        // 028
      " P 2c -2ac",      // 029
      " P 2 -2bc",       // 030
      " P 2ac -2",       // 031
      " P 2 -2ab",       // 032
      " P 2c -2n",       // 033
      " P 2 -2n",        // 034
      " C 2 -2",         // 035
      " C 2c -2",        // 036
      " C 2 -2c",        // 037
      " A 2 -2",         // 038
      " A 2 -2b",        // 039
      " A 2 -2a",        // 040
      " A 2 -2ab",       // 041
      " F 2 -2",         // 042
      " F 2 -2d",        // 043
      " I 2 -2",         // 044
      " I 2 -2c",        // 045
      " I 2 -2a",        // 046
      "-P 2 2",          // 047
      "-P 2ab 2bc",      // 048
      "-P 2 2c",         // 049
      "-P 2ab 2b",       // 050
      "-P 2a 2a",        // 051
      "-P 2a 2bc",       // 052
      "-P 2ac 2",        // 053
      "-P 2a 2ac",       // 054
      "-P 2 2ab",        // 055
      "-P 2ab 2ac",      // 056
      "-P 2c 2b",        // 057
      "-P 2 2n",         // 058
      "-P 2ab 2a",       // 059
      "-P 2n 2ab",       // 060
      "-P 2ac 2ab",      // 061
      "-P 2ac 2n",       // 062
      "-C 2c 2",         // 063
      "-C 2ac 2",        // 064
      "-C 2 2",          // 065
      "-C 2 2c",         // 066
      "-C 2a 2",         // 067
      "-C 2a 2ac",       // 068
      "-F 2 2",          // 069
      "-F 2uv 2vw",      // 070
      "-I 2 2",          // 071
      "-I 2 2c",         // 072
      "-I 2b 2c",        // 073
      "-I 2b 2",         // 074
      " P 4",            // 075
      " P 4w",           // 076
      " P 4c",           // 077
      " P 4cw",          // 078
      " I 4",            // 079
      " I 4bw",          // 080
      " P -4",           // 081
      " I -4",           // 082
      "-P 4",            // 083
      "-P 4c",           // 084
      "-P 4a",           // 085
      "-P 4bc",          // 086
      "-I 4",            // 087
      "-I 4ad",          // 088
      " P 4 2",          // 089
      " P 4ab 2ab",      // 090
      " P 4w 2c",        // 091
      " P 4abw 2nw",     // 092
      " P 4c 2",         // 093
      " P 4n 2n",        // 094
      " P 4cw 2c",       // 095
      " P 4nw 2abw",     // 096
      " I 4 2",          // 097
      " I 4bw 2bw",      // 098
      " P 4 -2",         // 099
      " P 4 -2ab",       // 100
      " P 4c -2c",       // 101
      " P 4n -2n",       // 102
      " P 4 -2c",        // 103
      " P 4 -2n",        // 104
      " P 4c -2",        // 105
      " P 4c -2ab",      // 106
      " I 4 -2",         // 107
      " I 4 -2c",        // 108
      " I 4bw -2",       // 109
      " I 4bw -2c",      // 110
      " P -4 2",         // 111
      " P -4 2c",        // 112
      " P -4 2ab",       // 113
      " P -4 2n",        // 114
      " P -4 -2",        // 115
      " P -4 -2c",       // 116
      " P -4 -2ab",      // 117
      " P -4 -2n",       // 118
      " I -4 -2",        // 119
      " I -4 -2c",       // 120
      " I -4 2",         // 121
      " I -4 2bw",       // 122
      "-P 4 2",          // 123
      "-P 4 2c",         // 124
      "-P 4a 2b",        // 125
      "-P 4a 2bc",       // 126
      "-P 4 2ab",        // 127
      "-P 4 2n",         // 128
      "-P 4a 2a",        // 129
      "-P 4a 2ac",       // 130
      "-P 4c 2",         // 131
      "-P 4c 2c",        // 132
      "-P 4ac 2b",       // 133
      "-P 4ac 2bc",      // 134
      "-P 4c 2ab",       // 135
      "-P 4n 2n",        // 136
      "-P 4ac 2a",       // 137
      "-P 4ac 2ac",      // 138
      "-I 4 2",          // 139
      "-I 4 2c",         // 140
      "-I 4bd 2",        // 141
      "-I 4bd 2c",       // 142
      " P 3",            // 143
      " P 31",           // 144
      " P 32",           // 145
      " R 3",            // 146
      "-P 3",            // 147
      "-R 3",            // 148
      " P 3 2",          // 149
      " P 3 2\"",        // 150
      " P 31 2 (0 0 4)", // 151
      " P 31 2\"",       // 152
      " P 32 2 (0 0 2)", // 153
      " P 32 2\"",       // 154
      " R 3 2\"",        // 155
      " P 3 -2\"",       // 156
      " P 3 -2",         // 157
      " P 3 -2\"c",      // 158
      " P 3 -2c",        // 159
      " R 3 -2\"",       // 160
      " R 3 -2\"c",      // 161
      "-P 3 2",          // 162
      "-P 3 2c",         // 163
      "-P 3 2\"",        // 164
      "-P 3 2\"c",       // 165
      "-R 3 2\"",        // 166
      "-R 3 2\"c",       // 167
      " P 6",            // 168
      " P 61",           // 169
      " P 65",           // 170
      " P 62",           // 171
      " P 64",           // 172
      " P 6c",           // 173
      " P -6",           // 174
      "-P 6",            // 175
      "-P 6c",           // 176
      " P 6 2",          // 177
      " P 61 2 (0 0 5)", // 178
      " P 65 2 (0 0 1)", // 179
      " P 62 2 (0 0 4)", // 180
      " P 64 2 (0 0 2)", // 181
      " P 6c 2c",        // 182
      " P 6 -2",         // 183
      " P 6 -2c",        // 184
      " P 6c -2",        // 185
      " P 6c -2c",       // 186
      " P -6 2",         // 187
      " P -6c 2",        // 188
      " P -6 -2",        // 189
      " P -6c -2c",      // 190
      "-P 6 2",          // 191
      "-P 6 2c",         // 192
      "-P 6c 2",         // 193
      "-P 6c 2c",        // 194
      " P 2 2 3",        // 195
      " F 2 2 3",        // 196
      " I 2 2 3",        // 197
      " P 2ac 2ab 3",    // 198
      " I 2b 2c 3",      // 199
      "-P 2 2 3",        // 200
      "-P 2ab 2bc 3",    // 201
      "-F 2 2 3",        // 202
      "-F 2uv 2vw 3",    // 203
      "-I 2 2 3",        // 204
      "-P 2ac 2ab 3",    // 205
      "-I 2b 2c 3",      // 206
      " P 4 2 3",        // 207
      " P 4n 2 3",       // 208
      " F 4 2 3",        // 209
      " F 4d 2 3",       // 210
      " I 4 2 3",        // 211
      " P 4acd 2ab 3",   // 212
      " P 4bd 2ab 3",    // 213
      " I 4bd 2c 3",     // 214
      " P -4 2 3",       // 215
      " F -4 2 3",       // 216
      " I -4 2 3",       // 217
      " P -4n 2 3",      // 218
      " F -4a 2 3",      // 219
      " I -4bd 2c 3",    // 220
      "-P 4 2 3",        // 221
      "-P 4a 2bc 3",     // 222
      "-P 4n 2 3",       // 223
      "-P 4bc 2bc 3",    // 224
      "-F 4 2 3",        // 225
      "-F 4a 2 3",       // 226
      "-F 4vw 2vw 3",    // 227
      "-F 4ud 2vw 3",    // 228
      "-I 4 2 3",        // 229
      "-I 4bd 2c 3",     // 230
  };
  return table[number];
}

const MatrixGroup &matrix_group_code_table(int number) {
  static const MatrixGroup *const table[] = {
      &matrix_group::unknown,    // 0
      &matrix_group::code_1,     // 1
      &matrix_group::code_1b,    // 2
      &matrix_group::code_2,     // 3
      &matrix_group::code_2,     // 4
      &matrix_group::code_2,     // 5
      &matrix_group::code_m,     // 6
      &matrix_group::code_m,     // 7
      &matrix_group::code_m,     // 8
      &matrix_group::code_m,     // 9
      &matrix_group::code_2_m,   // 10
      &matrix_group::code_2_m,   // 11
      &matrix_group::code_2_m,   // 12
      &matrix_group::code_2_m,   // 13
      &matrix_group::code_2_m,   // 14
      &matrix_group::code_2_m,   // 15
      &matrix_group::code_222,   // 16
      &matrix_group::code_222,   // 17
      &matrix_group::code_222,   // 18
      &matrix_group::code_222,   // 19
      &matrix_group::code_222,   // 20
      &matrix_group::code_222,   // 21
      &matrix_group::code_222,   // 22
      &matrix_group::code_222,   // 23
      &matrix_group::code_222,   // 24
      &matrix_group::code_mm2,   // 25
      &matrix_group::code_mm2,   // 26
      &matrix_group::code_mm2,   // 27
      &matrix_group::code_mm2,   // 28
      &matrix_group::code_mm2,   // 29
      &matrix_group::code_mm2,   // 30
      &matrix_group::code_mm2,   // 31
      &matrix_group::code_mm2,   // 32
      &matrix_group::code_mm2,   // 33
      &matrix_group::code_mm2,   // 34
      &matrix_group::code_mm2,   // 35
      &matrix_group::code_mm2,   // 36
      &matrix_group::code_mm2,   // 37
      &matrix_group::code_mm2,   // 38
      &matrix_group::code_mm2,   // 39
      &matrix_group::code_mm2,   // 40
      &matrix_group::code_mm2,   // 41
      &matrix_group::code_mm2,   // 42
      &matrix_group::code_mm2,   // 43
      &matrix_group::code_mm2,   // 44
      &matrix_group::code_mm2,   // 45
      &matrix_group::code_mm2,   // 46
      &matrix_group::code_mmm,   // 47
      &matrix_group::code_mmm,   // 48
      &matrix_group::code_mmm,   // 49
      &matrix_group::code_mmm,   // 50
      &matrix_group::code_mmm,   // 51
      &matrix_group::code_mmm,   // 52
      &matrix_group::code_mmm,   // 53
      &matrix_group::code_mmm,   // 54
      &matrix_group::code_mmm,   // 55
      &matrix_group::code_mmm,   // 56
      &matrix_group::code_mmm,   // 57
      &matrix_group::code_mmm,   // 58
      &matrix_group::code_mmm,   // 59
      &matrix_group::code_mmm,   // 60
      &matrix_group::code_mmm,   // 61
      &matrix_group::code_mmm,   // 62
      &matrix_group::code_mmm,   // 63
      &matrix_group::code_mmm,   // 64
      &matrix_group::code_mmm,   // 65
      &matrix_group::code_mmm,   // 66
      &matrix_group::code_mmm,   // 67
      &matrix_group::code_mmm,   // 68
      &matrix_group::code_mmm,   // 69
      &matrix_group::code_mmm,   // 70
      &matrix_group::code_mmm,   // 71
      &matrix_group::code_mmm,   // 72
      &matrix_group::code_mmm,   // 73
      &matrix_group::code_mmm,   // 74
      &matrix_group::code_4,     // 75
      &matrix_group::code_4,     // 76
      &matrix_group::code_4,     // 77
      &matrix_group::code_4,     // 78
      &matrix_group::code_4,     // 79
      &matrix_group::code_4,     // 80
      &matrix_group::code_4b,    // 81
      &matrix_group::code_4b,    // 82
      &matrix_group::code_4_m,   // 83
      &matrix_group::code_4_m,   // 84
      &matrix_group::code_4_m,   // 85
      &matrix_group::code_4_m,   // 86
      &matrix_group::code_4_m,   // 87
      &matrix_group::code_4_m,   // 88
      &matrix_group::code_422,   // 89
      &matrix_group::code_422,   // 90
      &matrix_group::code_422,   // 91
      &matrix_group::code_422,   // 92
      &matrix_group::code_422,   // 93
      &matrix_group::code_422,   // 94
      &matrix_group::code_422,   // 95
      &matrix_group::code_422,   // 96
      &matrix_group::code_422,   // 97
      &matrix_group::code_422,   // 98
      &matrix_group::code_4mm,   // 99
      &matrix_group::code_4mm,   // 100
      &matrix_group::code_4mm,   // 101
      &matrix_group::code_4mm,   // 102
      &matrix_group::code_4mm,   // 103
      &matrix_group::code_4mm,   // 104
      &matrix_group::code_4mm,   // 105
      &matrix_group::code_4mm,   // 106
      &matrix_group::code_4mm,   // 107
      &matrix_group::code_4mm,   // 108
      &matrix_group::code_4mm,   // 109
      &matrix_group::code_4mm,   // 110
      &matrix_group::code_4b2m,  // 111
      &matrix_group::code_4b2m,  // 112
      &matrix_group::code_4b2m,  // 113
      &matrix_group::code_4b2m,  // 114
      &matrix_group::code_4bm2,  // 115
      &matrix_group::code_4bm2,  // 116
      &matrix_group::code_4bm2,  // 117
      &matrix_group::code_4bm2,  // 118
      &matrix_group::code_4bm2,  // 119
      &matrix_group::code_4bm2,  // 120
      &matrix_group::code_4b2m,  // 121
      &matrix_group::code_4b2m,  // 122
      &matrix_group::code_4_mmm, // 123
      &matrix_group::code_4_mmm, // 124
      &matrix_group::code_4_mmm, // 125
      &matrix_group::code_4_mmm, // 126
      &matrix_group::code_4_mmm, // 127
      &matrix_group::code_4_mmm, // 128
      &matrix_group::code_4_mmm, // 129
      &matrix_group::code_4_mmm, // 130
      &matrix_group::code_4_mmm, // 131
      &matrix_group::code_4_mmm, // 132
      &matrix_group::code_4_mmm, // 133
      &matrix_group::code_4_mmm, // 134
      &matrix_group::code_4_mmm, // 135
      &matrix_group::code_4_mmm, // 136
      &matrix_group::code_4_mmm, // 137
      &matrix_group::code_4_mmm, // 138
      &matrix_group::code_4_mmm, // 139
      &matrix_group::code_4_mmm, // 140
      &matrix_group::code_4_mmm, // 141
      &matrix_group::code_4_mmm, // 142
      &matrix_group::code_3,     // 143
      &matrix_group::code_3,     // 144
      &matrix_group::code_3,     // 145
      &matrix_group::code_3,     // 146
      &matrix_group::code_3b,    // 147
      &matrix_group::code_3b,    // 148
      &matrix_group::code_312,   // 149
      &matrix_group::code_321,   // 150
      &matrix_group::code_312,   // 151
      &matrix_group::code_321,   // 152
      &matrix_group::code_312,   // 153
      &matrix_group::code_321,   // 154
      &matrix_group::code_32,    // 155
      &matrix_group::code_3m1,   // 156
      &matrix_group::code_31m,   // 157
      &matrix_group::code_3m1,   // 158
      &matrix_group::code_31m,   // 159
      &matrix_group::code_3m,    // 160
      &matrix_group::code_3m,    // 161
      &matrix_group::code_3b1m,  // 162
      &matrix_group::code_3b1m,  // 163
      &matrix_group::code_3bm1,  // 164
      &matrix_group::code_3bm1,  // 165
      &matrix_group::code_3bm,   // 166
      &matrix_group::code_3bm,   // 167
      &matrix_group::code_6,     // 168
      &matrix_group::code_6,     // 169
      &matrix_group::code_6,     // 170
      &matrix_group::code_6,     // 171
      &matrix_group::code_6,     // 172
      &matrix_group::code_6,     // 173
      &matrix_group::code_6b,    // 174
      &matrix_group::code_6_m,   // 175
      &matrix_group::code_6_m,   // 176
      &matrix_group::code_622,   // 177
      &matrix_group::code_622,   // 178
      &matrix_group::code_622,   // 179
      &matrix_group::code_622,   // 180
      &matrix_group::code_622,   // 181
      &matrix_group::code_622,   // 182
      &matrix_group::code_6mm,   // 183
      &matrix_group::code_6mm,   // 184
      &matrix_group::code_6mm,   // 185
      &matrix_group::code_6mm,   // 186
      &matrix_group::code_6bm2,  // 187
      &matrix_group::code_6bm2,  // 188
      &matrix_group::code_6b2m,  // 189
      &matrix_group::code_6b2m,  // 190
      &matrix_group::code_6_mmm, // 191
      &matrix_group::code_6_mmm, // 192
      &matrix_group::code_6_mmm, // 193
      &matrix_group::code_6_mmm, // 194
      &matrix_group::code_23,    // 195
      &matrix_group::code_23,    // 196
      &matrix_group::code_23,    // 197
      &matrix_group::code_23,    // 198
      &matrix_group::code_23,    // 199
      &matrix_group::code_m3b,   // 200
      &matrix_group::code_m3b,   // 201
      &matrix_group::code_m3b,   // 202
      &matrix_group::code_m3b,   // 203
      &matrix_group::code_m3b,   // 204
      &matrix_group::code_m3b,   // 205
      &matrix_group::code_m3b,   // 206
      &matrix_group::code_432,   // 207
      &matrix_group::code_432,   // 208
      &matrix_group::code_432,   // 209
      &matrix_group::code_432,   // 210
      &matrix_group::code_432,   // 211
      &matrix_group::code_432,   // 212
      &matrix_group::code_432,   // 213
      &matrix_group::code_432,   // 214
      &matrix_group::code_4b3m,  // 215
      &matrix_group::code_4b3m,  // 216
      &matrix_group::code_4b3m,  // 217
      &matrix_group::code_4b3m,  // 218
      &matrix_group::code_4b3m,  // 219
      &matrix_group::code_4b3m,  // 220
      &matrix_group::code_m3bm,  // 221
      &matrix_group::code_m3bm,  // 222
      &matrix_group::code_m3bm,  // 223
      &matrix_group::code_m3bm,  // 224
      &matrix_group::code_m3bm,  // 225
      &matrix_group::code_m3bm,  // 226
      &matrix_group::code_m3bm,  // 227
      &matrix_group::code_m3bm,  // 228
      &matrix_group::code_m3bm,  // 229
      &matrix_group::code_m3bm,  // 230
  };
  return *table[number];
}

const NormalizerGenerators &normalizer_table(int number) {
  static const NormalizerGenerators table[] = {
      {nullptr, nullptr},   // 0
      {"-1", nullptr},      // 1
      {nullptr, nullptr},   // 2
      {"-1", nullptr},      // 3
      {"-1", nullptr},      // 4
      {"-1", nullptr},      // 5
      {"-1", nullptr},      // 6
      {"-1", nullptr},      // 7
      {"-1", nullptr},      // 8
      {"-1", nullptr},      // 9
      {nullptr, nullptr},   // 10
      {nullptr, nullptr},   // 11
      {nullptr, nullptr},   // 12
      {nullptr, nullptr},   // 13
      {nullptr, nullptr},   // 14
      {nullptr, nullptr},   // 15
      {"-1", "3* -2'"},     // 16
      {"-1", "-2'w"},       // 17
      {"-1", "-2'"},        // 18
      {"-1", "3* -2'd"},    // 19
      {"-1", "-2'w"},       // 20
      {"-1", "-2'"},        // 21
      {"-1", "3* -2'"},     // 22
      {"-1", "3* -2'"},     // 23
      {"-1", "3* -2'd"},    // 24
      {"-1", "-2'"},        // 25
      {"-1", nullptr},      // 26
      {"-1", "-2'"},        // 27
      {"-1", nullptr},      // 28
      {"-1", nullptr},      // 29
      {"-1", nullptr},      // 30
      {"-1", nullptr},      // 31
      {"-1", "-2'"},        // 32
      {"-1", nullptr},      // 33
      {"-1", "-2'"},        // 34
      {"-1", "-2'"},        // 35
      {"-1", nullptr},      // 36
      {"-1", "-2'"},        // 37
      {"-1", nullptr},      // 38
      {"-1", nullptr},      // 39
      {"-1", nullptr},      // 40
      {"-1", nullptr},      // 41
      {"-1", "-2'"},        // 42
      {"-1uv", "-2'"},      // 43
      {"-1", "-2'"},        // 44
      {"-1", "-2'"},        // 45
      {"-1", nullptr},      // 46
      {nullptr, "3* -2'"},  // 47
      {nullptr, "3* -2'"},  // 48
      {nullptr, "-2'"},     // 49
      {nullptr, "-2'"},     // 50
      {nullptr, nullptr},   // 51
      {nullptr, nullptr},   // 52
      {nullptr, nullptr},   // 53
      {nullptr, nullptr},   // 54
      {nullptr, "-2'"},     // 55
      {nullptr, "-2'"},     // 56
      {nullptr, nullptr},   // 57
      {nullptr, "-2'"},     // 58
      {nullptr, "-2'"},     // 59
      {nullptr, nullptr},   // 60
      {nullptr, "3*"},      // 61
      {nullptr, nullptr},   // 62
      {nullptr, nullptr},   // 63
      {nullptr, nullptr},   // 64
      {nullptr, "-2'"},     // 65
      {nullptr, "-2'"},     // 66
      {nullptr, "-2'buv"},  // 67
      {nullptr, "-2'buv"},  // 68
      {nullptr, "3* -2'"},  // 69
      {nullptr, "3* -2'"},  // 70
      {nullptr, "3* -2'"},  // 71
      {nullptr, "-2'"},     // 72
      {nullptr, "3* -2'd"}, // 73
      {nullptr, "-2'bd"},   // 74
      {"-1", "-2'"},        // 75
      {nullptr, "2\""},     // 76
      {"-1", "-2'"},        // 77
      {nullptr, "2\""},     // 78
      {"-1", "-2'"},        // 79
      {"-1a", "2\""},       // 80
      {"-1", "-2'"},        // 81
      {"-1", "-2'"},        // 82
      {nullptr, "-2'"},     // 83
      {nullptr, "-2'"},     // 84
      {nullptr, "-2'"},     // 85
      {nullptr, "-2'"},     // 86
      {nullptr, "-2'"},     // 87
      {nullptr, "-2'd"},    // 88
      {"-1", nullptr},      // 89
      {"-1", nullptr},      // 90
      {nullptr, nullptr},   // 91
      {nullptr, nullptr},   // 92
      {"-1", nullptr},      // 93
      {"-1", nullptr},      // 94
      {nullptr, nullptr},   // 95
      {nullptr, nullptr},   // 96
      {"-1", nullptr},      // 97
      {"-1aw", nullptr},    // 98
      {"-1", nullptr},      // 99
      {"-1", nullptr},      // 100
      {"-1", nullptr},      // 101
      {"-1", nullptr},      // 102
      {"-1", nullptr},      // 103
      {"-1", nullptr},      // 104
      {"-1", nullptr},      // 105
      {"-1", nullptr},      // 106
      {"-1", nullptr},      // 107
      {"-1", nullptr},      // 108
      {"-1a", nullptr},     // 109
      {"-1a", nullptr},     // 110
      {"-1", nullptr},      // 111
      {"-1", nullptr},      // 112
      {"-1", nullptr},      // 113
      {"-1", nullptr},      // 114
      {"-1", nullptr},      // 115
      {"-1", nullptr},      // 116
      {"-1", nullptr},      // 117
      {"-1", nullptr},      // 118
      {"-1", nullptr},      // 119
      {"-1", nullptr},      // 120
      {"-1", nullptr},      // 121
      {"-1aw", nullptr},    // 122
      {nullptr, nullptr},   // 123
      {nullptr, nullptr},   // 124
      {nullptr, nullptr},   // 125
      {nullptr, nullptr},   // 126
      {nullptr, nullptr},   // 127
      {nullptr, nullptr},   // 128
      {nullptr, nullptr},   // 129
      {nullptr, nullptr},   // 130
      {nullptr, nullptr},   // 131
      {nullptr, nullptr},   // 132
      {nullptr, nullptr},   // 133
      {nullptr, nullptr},   // 134
      {nullptr, nullptr},   // 135
      {nullptr, nullptr},   // 136
      {nullptr, nullptr},   // 137
      {nullptr, nullptr},   // 138
      {nullptr, nullptr},   // 139
      {nullptr, nullptr},   // 140
      {nullptr, nullptr},   // 141
      {nullptr, nullptr},   // 142
      {"-1", "2 -2'"},      // 143
      {nullptr, "2 2\""},   // 144
      {nullptr, "2 2\""},   // 145
      {"-1", "-2\""},       // 146
      {nullptr, "2 -2'"},   // 147
      {nullptr, "-2\""},    // 148
      {"-1", "2"},          // 149
      {"-1", "2"},          // 150
      {nullptr, "2"},       // 151
      {nullptr, "2"},       // 152
      {nullptr, "2"},       // 153
      {nullptr, "2"},       // 154
      {"-1", nullptr},      // 155
      {"-1", "2"},          // 156
      {"-1", "2"},          // 157
      {"-1", "2"},          // 158
      {"-1", "2"},          // 159
      {"-1", nullptr},      // 160
      {"-1", nullptr},      // 161
      {nullptr, "2"},       // 162
      {nullptr, "2"},       // 163
      {nullptr, "2"},       // 164
      {nullptr, "2"},       // 165
      {nullptr, nullptr},   // 166
      {nullptr, nullptr},   // 167
      {"-1", "-2'"},        // 168
      {nullptr, "2\""},     // 169
      {nullptr, "2\""},     // 170
      {nullptr, "2\""},     // 171
      {nullptr, "2\""},     // 172
      {"-1", "-2'"},        // 173
      {"-1", "-2'"},        // 174
      {nullptr, "-2'"},     // 175
      {nullptr, "-2'"},     // 176
      {"-1", nullptr},      // 177
      {nullptr, nullptr},   // 178
      {nullptr, nullptr},   // 179
      {nullptr, nullptr},   // 180
      {nullptr, nullptr},   // 181
      {"-1", nullptr},      // 182
      {"-1", nullptr},      // 183
      {"-1", nullptr},      // 184
      {"-1", nullptr},      // 185
      {"-1", nullptr},      // 186
      {"-1", nullptr},      // 187
      {"-1", nullptr},      // 188
      {"-1", nullptr},      // 189
      {"-1", nullptr},      // 190
      {nullptr, nullptr},   // 191
      {nullptr, nullptr},   // 192
      {nullptr, nullptr},   // 193
      {nullptr, nullptr},   // 194
      {"-1", "-2'"},        // 195
      {"-1", "-2'"},        // 196
      {"-1", "-2'"},        // 197
      {"-1", "-2'd"},       // 198
      {"-1", "-2'd"},       // 199
      {nullptr, "-2'"},     // 200
      {nullptr, "-2'"},     // 201
      {nullptr, "-2'"},     // 202
      {nullptr, "-2'"},     // 203
      {nullptr, "-2'"},     // 204
      {nullptr, nullptr},   // 205
      {nullptr, "-2'd"},    // 206
      {"-1", nullptr},      // 207
      {"-1", nullptr},      // 208
      {"-1", nullptr},      // 209
      {"-1d", nullptr},     // 210
      {"-1", nullptr},      // 211
      {nullptr, nullptr},   // 212
      {nullptr, nullptr},   // 213
      {"-1", nullptr},      // 214
      {"-1", nullptr},      // 215
      {"-1", nullptr},      // 216
      {"-1", nullptr},      // 217
      {"-1", nullptr},      // 218
      {"-1", nullptr},      // 219
      {"-1", nullptr},      // 220
      {nullptr, nullptr},   // 221
      {nullptr, nullptr},   // 222
      {nullptr, nullptr},   // 223
      {nullptr, nullptr},   // 224
      {nullptr, nullptr},   // 225
      {nullptr, nullptr},   // 226
      {nullptr, nullptr},   // 227
      {nullptr, nullptr},   // 228
      {nullptr, nullptr},   // 229
      {nullptr, nullptr},   // 230
  };
  return table[number];
}

const std::vector<TwoFold> &reduced_cell_two_folds() {
  static const std::vector<TwoFold> table = {
      {{-1, -1, -1, 0, 0, 1, 0, 1, 0}, {-1, 1, 1}, {0, 1, 1}},
      {{-1, -1, 0, 0, 1, 0, 0, -1, -1}, {1, -2, 1}, {0, 1, 0}},
      {{-1, -1, 0, 0, 1, 0, 0, 0, -1}, {-1, 2, 0}, {0, 1, 0}},
      {{-1, -1, 0, 0, 1, 0, 0, 1, -1}, {-1, 2, 1}, {0, 1, 0}},
      {{-1, -1, 1, 0, 0, -1, 0, -1, 0}, {1, -1, 1}, {0, -1, 1}},
      {{-1, 0, -1, 0, -1, -1, 0, 0, 1}, {-1, -1, 2}, {0, 0, 1}},
      {{-1, 0, -1, 0, -1, 0, 0, 0, 1}, {-1, 0, 2}, {0, 0, 1}},
      {{-1, 0, -1, 0, -1, 1, 0, 0, 1}, {-1, 1, 2}, {0, 0, 1}},
      {{-1, 0, 0, -1, 0, -1, 1, -1, 0}, {0, -1, 1}, {1, -1, 1}},
      {{-1, 0, 0, -1, 0, 1, -1, 1, 0}, {0, 1, 1}, {-1, 1, 1}},
      {{-1, 0, 0, -1, 1, -1, 0, 0, -1}, {0, 1, 0}, {1, -2, 1}},
      {{-1, 0, 0, -1, 1, 0, 0, 0, -1}, {0, 1, 0}, {-1, 2, 0}},
      {{-1, 0, 0, -1, 1, 1, 0, 0, -1}, {0, 1, 0}, {-1, 2, 1}},
      {{-1, 0, 0, 0, -1, -1, 0, 0, 1}, {0, -1, 2}, {0, 0, 1}},
      {{-1, 0, 0, 0, -1, 0, -1, -1, 1}, {0, 0, 1}, {-1, -1, 2}},
      {{-1, 0, 0, 0, -1, 0, -1, 0, 1}, {0, 0, 1}, {-1, 0, 2}},
      {{-1, 0, 0, 0, -1, 0, -1, 1, 1}, {0, 0, 1}, {-1, 1, 2}},
      {{-1, 0, 0, 0, -1, 0, 0, -1, 1}, {0, 0, 1}, {0, -1, 2}},
      {{-1, 0, 0, 0, -1, 0, 0, 0, 1}, {0, 0, 1}, {0, 0, 1}},
      {{-1, 0, 0, 0, -1, 0, 0, 1, 1}, {0, 0, 1}, {0, 1, 2}},
      {{-1, 0, 0, 0, -1, 0, 1, -1, 1}, {0, 0, 1}, {1, -1, 2}},
      {{-1, 0, 0, 0, -1, 0, 1, 0, 1}, {0, 0, 1}, {1, 0, 2}},
      {{-1, 0, 0, 0, -1, 0, 1, 1, 1}, {0, 0, 1}, {1, 1, 2}},
      {{-1, 0, 0, 0, -1, 1, 0, 0, 1}, {0, 1, 2}, {0, 0, 1}},
      {{-1, 0, 0, 0, 0, -1, 0, -1, 0}, {0, -1, 1}, {0, -1, 1}},
      {{-1, 0, 0, 0, 0, 1, 0, 1, 0}, {0, 1, 1}, {0, 1, 1}},
      {{-1, 0, 0, 0, 1, -1, 0, 0, -1}, {0, 1, 0}, {0, -2, 1}},
      {{-1, 0, 0, 0, 1, 0, 0, -1, -1}, {0, -2, 1}, {0, 1, 0}},
      {{-1, 0, 0, 0, 1, 0, 0, 0, -1}, {0, 1, 0}, {0, 1, 0}},
      {{-1, 0, 0, 0, 1, 0, 0, 1, -1}, {0, 2, 1}, {0, 1, 0}},
      {{-1, 0, 0, 0, 1, 1, 0, 0, -1}, {0, 1, 0}, {0, 2, 1}},
      {{-1, 0, 0, 1, 0, -1, -1, -1, 0}, {0, -1, 1}, {-1, -1, 1}},
      {{-1, 0, 0, 1, 0, 1, 1, 1, 0}, {0, 1, 1}, {1, 1, 1}},
      {{-1, 0, 0, 1, 1, -1, 0, 0, -1}, {0, 1, 0}, {-1, -2, 1}},
      {{-1, 0, 0, 1, 1, 0, 0, 0, -1}, {0, 1, 0}, {1, 2, 0}},
      {{-1, 0, 0, 1, 1, 1, 0, 0, -1}, {0, 1, 0}, {1, 2, 1}},
      {{-1, 0, 1, 0, -1, -1, 0, 0, 1}, {1, -1, 2}, {0, 0, 1}},
      {{-1, 0, 1, 0, -1, 0, 0, 0, 1}, {1, 0, 2}, {0, 0, 1}},
      {{-1, 0, 1, 0, -1, 1, 0, 0, 1}, {1, 1, 2}, {0, 0, 1}},
      {{-1, 1, -1, 0, 0, -1, 0, -1, 0}, {-1, -1, 1}, {0, -1, 1}},
      {{-1, 1, 0, 0, 1, 0, 0, -1, -1}, {-1, -2, 1}, {0, 1, 0}},
      {{-1, 1, 0, 0, 1, 0, 0, 0, -1}, {1, 2, 0}, {0, 1, 0}},
      {{-1, 1, 0, 0, 1, 0, 0, 1, -1}, {1, 2, 1}, {0, 1, 0}},
      {{-1, 1, 1, 0, 0, 1, 0, 1, 0}, {1, 1, 1}, {0, 1, 1}},
      {{0, -1, -1, -1, 0, 1, 0, 0, -1}, {-1, 1, 0}, {-1, 1, 1}},
      {{0, -1, -1, 0, -1, 0, -1, 1, 0}, {-1, 0, 1}, {-1, 1, 1}},
      {{0, -1, 0, -1, 0, 0, -1, 1, -1}, {-1, 1, 1}, {-1, 1, 0}},
      {{0, -1, 0, -1, 0, 0, 0, 0, -1}, {-1, 1, 0}, {-1, 1, 0}},
      {{0, -1, 0, -1, 0, 0, 1, -1, -1}, {1, -1, 1}, {-1, 1, 0}},
      {{0, -1, 1, -1, 0, -1, 0, 0, -1}, {-1, 1, 0}, {1, -1, 1}},
      {{0, -1, 1, 0, -1, 0, 1, -1, 0}, {1, 0, 1}, {1, -1, 1}},
      {{0, 0, -1, -1, -1, 1, -1, 0, 0}, {-1, 1, 1}, {-1, 0, 1}},
      {{0, 0, -1, 0, -1, 0, -1, 0, 0}, {-1, 0, 1}, {-1, 0, 1}},
      {{0, 0, -1, 1, -1, -1, -1, 0, 0}, {-1, -1, 1}, {-1, 0, 1}},
      {{0, 0, 1, -1, -1, -1, 1, 0, 0}, {1, -1, 1}, {1, 0, 1}},
      {{0, 0, 1, 0, -1, 0, 1, 0, 0}, {1, 0, 1}, {1, 0, 1}},
      {{0, 0, 1, 1, -1, 1, 1, 0, 0}, {1, 1, 1}, {1, 0, 1}},
      {{0, 1, -1, 0, -1, 0, -1, -1, 0}, {-1, 0, 1}, {-1, -1, 1}},
      {{0, 1, -1, 1, 0, -1, 0, 0, -1}, {1, 1, 0}, {-1, -1, 1}},
      {{0, 1, 0, 1, 0, 0, -1, -1, -1}, {-1, -1, 1}, {1, 1, 0}},
      {{0, 1, 0, 1, 0, 0, 0, 0, -1}, {1, 1, 0}, {1, 1, 0}},
      {{0, 1, 0, 1, 0, 0, 1, 1, -1}, {1, 1, 1}, {1, 1, 0}},
      {{0, 1, 1, 0, -1, 0, 1, 1, 0}, {1, 0, 1}, {1, 1, 1}},
      {{0, 1, 1, 1, 0, 1, 0, 0, -1}, {1, 1, 0}, {1, 1, 1}},
      {{1, -1, -1, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {-2, 1, 1}},
      {{1, -1, 0, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {-2, 1, 0}},
      {{1, -1, 1, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {2, -1, 1}},
      {{1, 0, -1, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {-2, 0, 1}},
      {{1, 0, 0, -1, -1, 0, -1, 0, -1}, {-2, 1, 1}, {1, 0, 0}},
      {{1, 0, 0, -1, -1, 0, 0, 0, -1}, {-2, 1, 0}, {1, 0, 0}},
      {{1, 0, 0, -1, -1, 0, 1, 0, -1}, {2, -1, 1}, {1, 0, 0}},
      {{1, 0, 0, 0, -1, 0, -1, 0, -1}, {-2, 0, 1}, {1, 0, 0}},
      {{1, 0, 0, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {1, 0, 0}},
      {{1, 0, 0, 0, -1, 0, 1, 0, -1}, {2, 0, 1}, {1, 0, 0}},
      {{1, 0, 0, 1, -1, 0, -1, 0, -1}, {-2, -1, 1}, {1, 0, 0}},
      {{1, 0, 0, 1, -1, 0, 0, 0, -1}, {2, 1, 0}, {1, 0, 0}},
      {{1, 0, 0, 1, -1, 0, 1, 0, -1}, {2, 1, 1}, {1, 0, 0}},
      {{1, 0, 1, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {2, 0, 1}},
      {{1, 1, -1, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {-2, -1, 1}},
      {{1, 1, 0, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {2, 1, 0}},
      {{1, 1, 1, 0, -1, 0, 0, 0, -1}, {1, 0, 0}, {2, 1, 1}},
  };
  return table;
}

} // namespace mxi::stills::sgtbx
