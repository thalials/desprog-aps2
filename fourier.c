#include <math.h>

#include "fourier.h"

// nft -> implementa a formula normal
void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign)
{
    for (int k = 0; k < n; k++)
    {
        t[k] = 0;

        for (int j = 0; j < n; j++)
        {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n)
{
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n)
{
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++)
    {
        s[k] /= n;
    }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign)
{
    // declaracoes
    int metade = n / 2;
    double complex sp[metade];
    double complex si[metade];

    // tp e ti transformadas de Fourier de sp e si, respectivamente
    double complex tp[metade];
    double complex ti[metade];

    // parametro menor possivel -> caso em que o vetor é de tamanho unitario
    if (n == 1){
        t[0] = s[0];
        return;
    }

    for (int i = 0; i < metade; i++){
        sp[i] = s[i * 2];
        si[i] = s[i * 2 + 1];
    }

    // chamada recursiva, menor parametro -> metade
    fft(sp, tp, metade, sign);
    fft(si, ti, metade, sign);

    // para todo k de 0 a metade-1 -> t[k] = tp[k] + ti[k]*e^(2.pi.k.i/n)
    for (int k = 0; k < metade; k++){
        t[k] = tp[k] + ti[k] * cexp(sign * 2 * PI * k * I / n);
        t[k + metade] = tp[k] - ti[k] * cexp(sign * 2 * PI * k * I / n);
    }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n){
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n){
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++)
    {
        s[k] /= n;
    }
}

// height -> linhas
// width -> colunas

// como as ideias se repetem em fft_forward_2d() e fft_inverse_2d(), a função aplica_fft() abaixo foi criada
// para deixar o codigo menos poluido 
void aplica_fft(void fft_forward_inverse(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n),
                double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    
    double complex linha_s[width];
    double complex linha_t[width];

    for (int l = 0; l < height; l++){
        // colocando novos valores no vetor linha_s com elementos de matrix
        for (int c = 0; c < width; c++){
            linha_s[c] = matrix[l][c];
        }

        // transformada normal em linha_s
        // usei o valor de width pq o ultimo for assume no max, o valor witdh-1, porem ainda n está 100% claro
        fft_forward_inverse(linha_s, linha_t, width);

        // cada linha da matriz terá um valor da transformada 
        for (int c = 0; c < width; c++){
            matrix[l][c] = linha_t[c];
        }
    }

    // mesma ideia para o caso anterior ao aplicar mudancas nas colunas //
    double complex coluna_s[height];
    double complex coluna_t[height];

    for (int coluna = 0; coluna < width; coluna++) {
        // coloca novos valores no vetor coluna_s com elementos de matrix
        for (int linha = 0; linha < height; linha++) {
            coluna_s[linha] = matrix[linha][coluna];
        }

        // transformada normal em coluna_s
        fft_forward_inverse(coluna_s, coluna_t, height);

        // cada coluna da matriz terá um valor da transformada 
        for (int linha = 0; linha < height; linha++) {
            matrix[linha][coluna] = coluna_t[linha];
        }
    }
    return;
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height){   
    aplica_fft(fft_forward, matrix, width, height);
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height){
    aplica_fft(fft_inverse, matrix, width, height);
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip)
{
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip)
            {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height)
{
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height)
{
    filter(input, output, width, height, 1);
}
