/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1(){
	const double pkt = 62.744;
	matrix x(pkt);
	std::cout << "Wartosc funkcji celu dla punktu " << pkt << " = " << m2d(ff1T(x)) << "\n"; // wartosc minimum
	srand(time(0));
	double x0 = 0, alpha = 1.25, d = 2.0;
	int Nmax = 1000;
	double** p_wsk = new double*[100];
	std::cout << "\t=========================FIBONACI================================\n";
//	for (int i = 0; i < 100; i++) {
		//double x0 = rand() % 201 - 100; // losuje z przedzia³u [-100, 100]
	//	double* p = expansion(ff1T, x0, d, alpha, Nmax);		
		//p_wsk[i] = p;
		//std::cout << "," << x0 << "," << p[0] << "," << p[1] << "\n";
	//}
	std::cout << "\n" << p_wsk[0][0] << "," << p_wsk[0][1] << "\n";
	double a, b, epsilon = 0.1;
	//for (int i = 0; i < 100; i++) {
		//a = p_wsk[i][0];
		//b = p_wsk[i][1];
		if (a > b) {
			std::swap(a, b);
		}
		std::cout << a << "," << b << ",";
		solution xopt = fib(ff1T, -100, 100, epsilon);
		cout << m2d(xopt.x) << "," << m2d(xopt.y) << "\n";
	//}
	std::cout << "\t=======================LAGRANGE=============================\n";
	
	double c, dd, gamma = 1e-200;
	//for (int i = 0; i < 100; i++) {
	//	c = p_wsk[i][0];
		//dd = p_wsk[i][1];
	//	if (c > dd) {
	//		std::swap(c, dd);
		//}
	// 	solution xopt1 = lag(ff1T, -100, 100, epsilon, gamma, Nmax);
	 	//cout << m2d(xopt1.x) << "," << m2d(xopt1.y) << "\n";
	//}


}
void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
