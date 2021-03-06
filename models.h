//  models.h
//  xcssecovid
//
//  Created by Dr. Rolf Jansen on 2020-03-26.
//  Copyright © 2020 Dr. Rolf Jansen. All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
//  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
//  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
//  OF THE POSSIBILITY OF SUCH DAMAGE.


extern char *modelDescription_LF;
int initialValues_LF(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_LF(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_LDE;
int initialValues_LDE(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_LDE(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_SI;
int initialValues_SI(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_SI(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_SIR;
int initialValues_SIR(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_SIR(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_SEIR;
int initialValues_SEIR(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_SEIR(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_SEIR_de;
int initialValues_SEIR_de(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_SEIR_de(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_SIRX;
int initialValues_SIRX(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_SIRX(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_ERF;
int initialValues_ERF(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_ERF(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);

extern char *modelDescription_GLF;
int initialValues_GLF(ldouble t1, ldouble min, ldouble max, ldouble A[mpar], int f[mpar]);
int modelFunction_GLF(ldouble t, ldouble *Y, ldouble A[mpar], int f[mpar], bool fit, bool init);
