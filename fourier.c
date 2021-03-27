#include <math.h>

#include "fourier.h"

// nft -> implementa a formula normal
void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k] = 0;

        for (int j = 0; j < n; j++) {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    // declaracoes 
    int metade = n/2;
    double complex sp[ metade ];
    double complex si[ metade ];

    // tp e ti transformadas de Fourier de sp e si, respectivamente
    double complex tp[ metade ];
    double complex ti[ metade ];
    
    // parametro menor possivel -> caso em que o vetor é de tamanho unitario
    if (n == 1) {
        t[0] = s[0];
        return;
    }

    for (int i = 0; i < metade ; i++) {
        sp[i] = s[i*2];
        si[i] = s[i*2 + 1];
    }
   
    // chamada recursiva, menor parametro -> metade
    fft(sp, tp, metade, sign);
    fft(si, ti, metade, sign);

    // para todo k de 0 a metade - 1 -> t[k] = tp[k] + ti[k]*e^(2.pi.k.i/n)
    for (int k = 0; k < metade; k++ ) {
        t[k] = tp[k] + ti[k]*cexp(sign * 2 * PI * k * I / n);
        t[k + metade] = tp[k] - ti[k]*cexp(sign * 2 * PI * k * I / n);
    }

}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

//double complex matrix[metade][metade];

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    for(int l = 0; l < MAX_SIZE; l++){
        fft_forward(matrix[l], matrix, height);
    }
    for(int c = 0; c < MAX_SIZE; c++){
        fft_forward(matrix[c], matrix, width);
    }
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    for(int l = 0; l < MAX_SIZE; l++){
        fft_inverse(matrix[l], matrix, height);
    }
    for(int c = 0; c < MAX_SIZE; c++){
        fft_inverse(matrix[c], matrix, width);
    }
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
