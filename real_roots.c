#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define h 1e-6
#define Abs(x) ((x)>=0?(x):-(x))
static int total_soln;

static double f(const double x) {
    // write your function here!!
    return (x*x*x*x+0.5*x*x*x+5*x*x+15*x+3.587);
}

static double df(const double x) {
    return ((f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))/(12*h));
}

static double ddf(const double x) {
    return ((-f(x-2*h)+16*f(x-h)-30*f(x)+16*f(x+h)-f(x+2*h))/(12*h*h));
}

static double newton_solve(const double x0) {
    if (Abs(df(x0))>1e-9 && Abs(ddf(x0)/(2*df(x0)))<1) {
        double err=100;
        double xi=x0,xi_1=0;
        int it=0;
        while (err>1e-8) {
            if (it==100000 && err>1e-8) {
                printf("Didn't converge! Try another guess.\n");
                return NAN;
            }
            if (Abs(df(xi))<1e-9) {
                printf("Cannot solve!\nf'(x)~0! Try another guess.\n");
                return NAN;
            }
            xi_1=xi-f(xi)/df(xi);
            err=Abs(xi-xi_1);
            xi=xi_1;
            it++;
        }
        return xi_1;
    }
    printf("Cannot solve! Try another guess.\n");
    return NAN;
}

static void real_solve(const double l, const double u,double **soln) {
    double *x=malloc(1e6*sizeof(double));
    if (!x) {
        return;
    }
    double *y=malloc(1e6*sizeof(double));
    if (!y) {
        free(x);
        return;
    }

    double d=(u-l)/(1e6-1);
    for (int i=0; i<1e6; i++) {
        x[i]=l+i*d;
        y[i]=f(x[i]);
    }

    int no_soln=1;
    for (int i=1; i<1e6; i++) {
        if (y[0]*y[i]<0) {
            no_soln=0;
            break;
        }
    }
    if (no_soln) {
        printf("No real solution exists in this range!\n");
        free(x);
        free(y);
        return;
    }

    int sz=0;
    double low= newton_solve(l);
    if (isnan(low)) {
        printf("Please change the lower range! Not converging.\n");
        free(x);
        free(y);
        return;
    }
    double up= newton_solve(u);
    if (isnan(up)) {
        printf("Please change the upper range! Not converging.\n");
        free(x);
        free(y);
        return;
    }
    if (Abs(low-up)<1e-8) {
        *soln=malloc(sizeof(double));
        total_soln=1;
        if (!(*soln)) {
            free(x);
            free(y);
        }
        (*soln)[0]=low;
        return;
    }
    d=(up-low)/(1e6-1);
    for (int i=0; i<1e6; i++) {
        x[i]=low+i*d;
        y[i]=f(x[i]);
    }
    for (int i=0; i<1e6; i++) {
        if (Abs(y[i])<1e-8) {
            sz++;
        }
    }
    total_soln=sz;
    *soln=malloc(sizeof(double)*sz);
    if (!soln) {
        free(x);
        free(y);
    }
    int it=0;
    for (int i=0; i<1e6; i++) {
        if (Abs(y[i])<1e-8) {
            (*soln)[it++]=x[i];
        }
    }
    free(x);
    free(y);
}

int main() {
    double *soln=NULL;
    real_solve(-1050.254,25.541,&soln);

    // printf("%d\n",total_soln);
    for (int i=0; i<total_soln; i++) {
        printf("%.10lf\n",soln[i]);
    }

    free(soln);
}