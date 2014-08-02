#ifndef QUADFILTER_H
#define QUADFILTER_H
#include <math.h>
#include <assert.h>

// A quadratic IIR discrete time filter
class QuadFilter {
    double c_0;
    double c_1;
    double c_2;
    double d_1;
    double d_2;
    double input_1;
    double input_2;
    double output_1;
    double output_2;

public:
    static QuadFilter CreateBandpass(double f_min, double f_max,
            double f_sampling, double input_start) {
        QuadFilter filter;

        double a = tan(M_PI * f_min * (1.0/f_sampling));
        double b = tan(M_PI * f_max * (1.0/f_sampling));
        filter.c_0 = -b/((1+a)*(1+b));
        filter.c_1 = 0;
        filter.c_2 = -filter.c_0;
        filter.d_1 = (2 - 2*a*b)/((1+a)*(1+b));
        filter.d_2 = -((1-a)*(1-b))/((1+a)*(1+b));
        filter.input_1 = filter.input_2 = input_start;
        filter.output_1 = filter.output_2 = 0;

        return filter;
    }

    static QuadFilter CreateHighpass(double f_cutoff, unsigned int order,
            unsigned int iteration, double f_sampling, double input_start) {
        QuadFilter filter;

	// we only support even orders at this time
	// the way to support odd order filters might be chaining
	// quadratic filters and a linear filter
	assert(order % 2 == 0);

        double a = tan(M_PI * f_cutoff * (1.0/f_sampling));
	double alpha = (a*a);
	double gamma_j = 2 * a * sin(M_PI * (2+(2*iteration))/(2*order));
	double denominator = 1 + alpha + gamma_j;

        filter.c_0 = 1/denominator;
        filter.c_1 = -2/denominator;
        filter.c_2 = 1/denominator;
        filter.c_1 = (2 - (2*alpha))/denominator;
        filter.c_2 = (gamma_j - alpha - 1)/denominator;
        filter.input_1 = filter.input_2 = input_start;
        filter.output_1 = filter.output_2 = 0;

        return filter;
    }

    static QuadFilter CreateLowpass(double f_cutoff, unsigned int order,
            unsigned int iteration, double f_sampling, double input_start) {
        QuadFilter filter;

	// we only support even orders at this time
	// the way to support odd order filters might be chaining
	// quadratic filters and a linear filter
	assert(order % 2 == 0);

        double a = tan(M_PI * f_cutoff * (1.0/f_sampling));
	double alpha = (a*a);
	double gamma_j = 2 * a * sin(M_PI * (2+(2*iteration))/(2*order));
	double denominator = 1 + alpha + gamma_j;

        filter.c_0 = -alpha/denominator;
        filter.c_1 = -(2*alpha)/denominator;
        filter.c_2 = -alpha/denominator;
        filter.c_1 = (2 - (2*alpha))/denominator;
        filter.c_2 = (gamma_j - alpha - 1)/denominator;
        filter.input_1 = filter.input_2 = input_start;
        filter.output_1 = filter.output_2 = input_start;

        return filter;
    }

    double filterInput(double input_0) {
        // calculate the new output
        double output_0 = c_0 * input_0 + c_1 * input_1 + c_2 * input_2 +
            d_1 * output_1 + d_2 * output_2;

        // shift the new data into the previous values
        input_2 = input_1;
        input_1 = input_0;
        output_2 = output_1;
        output_1 = output_0;

        return output_0;
    }
};

#endif // QUADFILTER_H
