#include <math.h>
#include <cstdio>

/* function declaration - network evaluation */
void export_eval(const double*, double*);


/* ------------------------------------------------------------ */


static const double export_in_A[9] = {
 0.005871191574553051, 0, 0, 0,
 0.890015910898966, 0, 0, 0,
 0.4345734771979852
};


static const double export_in_b[3] = {
 -2.530859438885565, -1.011336515513126, -0.9054317306724021
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
 4.27192764738202
};


static const double export_out_b[1] = {
 3.655489332073003
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
 0.7694799029706706, 0.0009323679992975241, -0.009673029928858587, 1.142485804830319,
 0.2596257753819519, 0.5703600579613278, 0.3216019626510759, 0.3412842062772985,
 0.288225347618272, -0.2156982791365307, -0.2739939219324733, 0.119991776776743,
 0.3097299967059008, 0.2327818961864505, 0.2165917711199529, 0.01948695099226739,
 0.1290976921008944, -1.748262681700963, -0.4494221576904884, 0.07905987847189418,
 0.001035772595245012, -0.1819797605028023, 0.4763010121474991, -0.01629637924856567,
 0.332088409758003, 0.1919871015634123, -0.2934852336152938, 0.4897192633227027,
 -0.1973425733648359, 0.2982758693594617, 0.2378391666236043, 0.2246231687638703,
 0.1268111461135582, 0.3323204636502285, 0.3609274316533565, 1.171258375919519,
 0.1994576288519849, 0.6474501465455454, 0.05706473382971627, 0.3646895505619482,
 0.1075026920937489, -0.6976416266778649, 0.05561742189412289, 0.6021062463973895,
 4.894655113647565, -0.2395449305687137, -0.2031485693875097, 0.4844236482255392,
 0.2423964820958496, -0.4701803826034681, 0.5042865281748397, 0.1155220150638223,
 -0.5389050098928198, -0.2394297108458284, -0.1844930089152931, -0.667424714393407,
 -0.2780162586514213, -1.553209024771389, -0.9999250999470504, -0.4021890588283863,
 2.000660652590356
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
    Below is a test for H2-air, stoichiometric, 1 bar and 300 K
     */
    double testinp[3] = {300, 1.0, 1.0};
    double testout = 0.0;

    double *inp = testinp;
    double *out = &testout;

    export_eval(inp, out);

    printf("LBV is %f m/s\n", *out);

    return 0;
}
