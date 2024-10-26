#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix funkcja_testowa_lab1(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	double PI = 3.14;
	y = -cos(0.1 * m2d(x)) * exp(-pow(0.1 * m2d(x) - 2 * PI, 2)) + 0.002 * pow(0.1 * x, 2);
	return y;
}

double GetFib(int n) {
	if (n <= 1)
		return n;

	double prev1 = 0;
	double prev2 = 1;
	double current = 0;

	for (int i = 2; i <= n; ++i) {
		current = prev1 + prev2;
		prev1 = prev2;
		prev2 = current;
	}

	return current;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(3, 1);
	double a = 0.98, b = 0.63, g = 9.81;

	double PA = 0.5;
	double TA = 90.0;

	double PB = 1.0;
	double DB = 0.00365665;

	double Fin = 0.01;
	double Tin = 20.0;

	double Aout = Y(0) > 0 ? a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA) : 0;
	double Bout = Y(1) > 1 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;

	dY(0) = -Aout;
	dY(1) = Aout + Fin - Bout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + Aout / Y(1) * (TA - Y(2));

	return dY;
}

matrix fun_rzeczywista_lab1(matrix x, matrix ud1, matrix ud2)
{
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});

	double t0 = 0.0;
	double tend = 2000.0;
	double dt = 1.0;

	matrix* Y = solve_ode(df1, t0, dt, tend, Y0, ud1, x);
	double max = Y[1](0, 2);

	int length = get_len(Y[0]);
	double y;
	for (int i = 1; i < length; i++)
	{
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	}
	y = abs(max - 50);
	return y;
}