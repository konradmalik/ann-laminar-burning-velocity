#include <math.h>
#include <cstdio>

/* function declaration - network evaluation */
void export_eval(const double*, double*);


/* ------------------------------------------------------------ */


static const double export_in_A[9] = {
 0.003247474311062158, 0, 0, 0,
 0.08542330830352415, 0, 0, 0,
 1.656291206135793
};


static const double export_in_b[3] = {
 -1.801080678029856, -0.8752447691403317, -1.817307517333532
};


void
export_in(const double *x, double *y)
{
    int k = 0, i, j;
    for (i = 0; i < 3; ++i) {
        y[i] = export_in_b[i];
        for (j = 0; j < 3; ++j, ++k) {
            y[i] += export_in_A[k] * x[j];
        }
    }
}


static const double export_out_A[1] = {
 1.484687479545055
};


static const double export_out_b[1] = {
 1.261964407955686
};


void
export_out(const double *x, double *y)
{
    int k = 0, i, j;
    for (i = 0; i < 1; ++i) {
        y[i] = export_out_b[i];
        for (j = 0; j < 1; ++j, ++k) {
            y[i] += export_out_A[k] * x[j];
        }
    }
}


static const double export_w[61] = {
 -0.325179904501508, 0.3630731653787407, -0.3171570519909662, -0.8896395318625778,
 -0.1964710751553572, 0.2720470390794683, 0.1749395067111781, -0.2037283290077966,
 0.09787314703797428, -0.9333946214951209, 1.246035782993538, -0.9803404839873349,
 -0.1963679337187418, -0.5107281159998256, 1.521480827630051, -2.326068933127351,
 0.06213266802348242, -0.2770552520323777, 0.6073577593754123, -0.2513922582275077,
 -0.0457368471529099, -0.08807329377085521, 0.0392923747060681, -0.7739184469122425,
 0.3120011532098662, 0.2099986433006529, 0.6194663344236687, -1.681655661338133,
 -0.01145554040443144, 0.6820504116306647, 1.014711920263665, 0.1560157964269563,
 -0.536919938153857, -0.2203001640521614, 1.16514649662503, 0.6621162348032152,
 2.575519009122875, -8.820743852808572, 5.23200843980797, 3.643148631144902,
 -4.612957904056657, -0.1646073412170468, -0.03610102760127006, 0.647141245440874,
 0.4608211696472901, -0.1725339678438308, -0.1980742704545649, 0.03058695440954697,
 0.003980053839517748, 0.07147530839732344, -0.09486893321070347, -0.02136187713576113,
 -0.8104499725509378, -0.3475306519232343, -0.04453414375872867, -0.2822147292125957,
 -0.4200143781739888, -0.363548685471267, -1.435533016652866, -2.13810709537974,
 -0.7563985027801263
};


void
export_eval(const double *input, double *output)
{
    double n[16];
    double x;
    int k = 0, i = 3, j;

    export_in(input, n);
    for (; i < 7; ++i) {
        x = export_w[k++];
        for (j = 0; j < 3; ++j) x += n[j] * export_w[k++];
        n[i] = (double)2. / ((double)1. + exp((double)-2 * x)) - (double)1.;
    }
    for (; i < 11; ++i) {
        x = export_w[k++];
        for (j = 3; j < 7; ++j) x += n[j] * export_w[k++];
        n[i] = (double)2. / ((double)1. + exp((double)-2 * x)) - (double)1.;
    }
    for (; i < 15; ++i) {
        x = export_w[k++];
        for (j = 7; j < 11; ++j) x += n[j] * export_w[k++];
        n[i] = (double)2. / ((double)1. + exp((double)-2 * x)) - (double)1.;
    }
    for (; i < 16; ++i) {
        x = export_w[k++];
        for (j = 11; j < 15; ++j) x += n[j] * export_w[k++];
        n[i] = x;
    }
    export_out(n + 15, output);
}


int main()
{
    /* inputs order:
    temperature [K], pressure [bar], equivalence ratio (phi)
    output:
    laminar burning velocity [m/s]
    Below is a test for C3H8-air, stoichiometric, 1 bar and 300 K
     */
    double testinp[3] = {300, 1.0, 1.0};
    double testout = 0.0;

    double *inp = testinp;
    double *out = &testout;

    export_eval(inp, out);

    printf("LBV is %f m/s\n", *out);

    return 0;
}
