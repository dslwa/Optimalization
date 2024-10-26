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
	double epsilon = 0.01;
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

void lab1()
{
	srand(time(NULL));

	double* ekspansja = new double[2];
	int Nmax = 1000;
	double epsilon = 0.001, gamma = 0.1;

	// Dane to tabeli 1
	ofstream ekspansja_tabela_1("./dane/lab_01/funkcja_testowa/exp_tab_1.txt");
	ofstream fibonaci_tabela_1("./dane/lab_01/funkcja_testowa/fib_tab_1.txt");
	ofstream lagrange_tabela_1("./dane/lab_01/funkcja_testowa/lag_tab_1.txt");

	// Trzy wpó³czynniki alfa dla ekspansji
	double alpha, alpha_1 = 1.3, alpha_2 = 1.6, alpha_3 = 2.3;
	double x;
	double d = 2.0;
	alpha = alpha_1;

	for (int i = 0; i < 300; i++)
	{
		int range = 800;  // Liczba mo¿liwych wartoœci (od -100 do 100 z krokiem 0.25 to 800 wartoœci)
		double x = -100 + (rand() % range) * 0.25;

		if (i == 100) alpha = alpha_2;
		else if (i == 200) alpha = alpha_3;
		ekspansja = expansion(funkcja_testowa_lab1, x, d, alpha, Nmax);
		ekspansja_tabela_1 << x << "," << ekspansja[0] << "," << ekspansja[1] << "," << solution::f_calls << "\n";
		std::cout <<  i  << " " << x << "," << ekspansja[0] << "," << ekspansja[1] << "," << solution::f_calls << "\n";
		solution::clear_calls();

		solution fibonaci_1 = fib(funkcja_testowa_lab1, ekspansja[0], ekspansja[1], epsilon);
		fibonaci_tabela_1 << m2d(fibonaci_1.x) << "," << m2d(fibonaci_1.y) << "," << solution::f_calls << "," << fibonaci_1.flag << "\n";
		solution::clear_calls();

		solution lagrange_1 = lag(funkcja_testowa_lab1, ekspansja[0], ekspansja[1], epsilon, gamma, Nmax);
		lagrange_tabela_1 << m2d(lagrange_1.x) << "," << m2d(lagrange_1.y) << "," << solution::f_calls << "," << lagrange_1.flag << "\n";
		solution::clear_calls();
	}
	ekspansja_tabela_1.close();
	fibonaci_tabela_1.close();
	lagrange_tabela_1.close();

	// Dane do wykresu
	ofstream fib_wykres("./dane/lab_01/funkcja_testowa/fibonaci_wykres.txt");
	ofstream lag_wykres("./dane/lab_01/funkcja_testowa/lagrange_wykres.txt");
	
	solution fib2 = fib(funkcja_testowa_lab1, -100, 100, epsilon);
	fib_wykres << fib2 << "\n\n" << fib2.ud << "\n";
	solution::clear_calls();

	solution lag2 = lag(funkcja_testowa_lab1, -100, 100, epsilon, gamma, Nmax);
	lag_wykres << lag2 << "\n\n" << lag2.ud << "\n";
	solution::clear_calls();

	fib_wykres.close();
	lag_wykres.close();

	// Dane do tabeli 3
	ofstream fib_tab_3("./dane/lab_01/funkcja_rzeczywista/fibonaci_tab_3.txt");
	ofstream lag_tab_3("./dane/lab_01/funkcja_rzeczywista/lagrange_tab_3.txt");

	solution fib3 = fib(fun_rzeczywista_lab1, 0.0001, 0.01, epsilon);
	fib_tab_3 << fib3 << endl;
	solution::clear_calls();

	solution lag3 = lag(fun_rzeczywista_lab1, 0.0001, 0.01, epsilon, gamma, Nmax);
	lag_tab_3 << lag3 << endl;
	solution::clear_calls();

	fib_tab_3.close();
	lag_tab_3.close();

	// Dane do symulacji
	ofstream fib_sym("./dane/lab_01/funkcja_rzeczywista/fibonaci_symulacja.csv");
	ofstream lag_sym("./dane/lab_01/funkcja_rzeczywista/lagrange_symulacja.csv");

	matrix Y0 = matrix(3, new double[3] {5, 1, 20});
	matrix* Y_fib = solve_ode(df1, 0, 1, 2000, Y0, NAN, fib3.x(0));
	solution::clear_calls();
	matrix* Y_lag = solve_ode(df1, 0, 1, 2000, Y0, NAN, lag3.x(0));

	fib_sym << Y_fib[1];
	lag_sym << Y_lag[1];

	fib_sym.close();
	lag_sym.close();
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
