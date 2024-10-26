#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}


double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{//funkcja celu, pocz�tkowa waro�c zmiennej, krok kt�ry okre�la odlweg�o�c mi�dzy punktem pocz�tkowaym a ko�cowym, alpha rozszerza lub zmienjsza krok z ka�dym kroku p�tli, max liczba wywo�a� funkcji celu
	try
	{
		double* p = new double[2] { 0, 0 }; //przechowuje nasze przedzia�y gdzie jest miinmum
		int i = 0;
		double temp = 0;
		solution X0(x0), X1(x0 + d); //punkt pocz�tkowy, obiekty klasy solution x1 = x0 + d
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);// wylicza warto�� funkcji celu w tych punktach
		if (X0.y == X1.y)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x); // Je�li te same wartosci to zwraca ten przedzia� bo zak�ada, �e nie ma w nim znaczacych zmien finkcji celu
			return p;
		}
		if (X1.y > X0.y)
		{
			d = -d; // je�li warto�c funcji w x1 > x0 to funkcja ro�nie czyli nie zmierzamy w stron� minimum, zmieniamy kierunek daltego d na - bo chcemy si� cofn��
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2); // liczymy funkcje celu 
			if (X1.y >= X0.y)
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x) - d; // je�li nadal jest wi�ksza to zwracamy ten przedzia
				return p;
			}
		}
		do
		{
			if (solution::f_calls > Nmax) //sprawdza czy liczba wywolan funkcji celu orzekroczy�a dozwoln� max liczb�
			{
				X0.flag = 0;
				return 0;
			}
			i = i + 1;
			temp = m2d(X0.x);

			X1 = X0.x + pow(alpha, i) * d; //oblicza nowy punkt x1 zwiekszajac krok alfa do poteki i
			X1.fit_fun(ff, ud1, ud2); // liczbymy wartosc funkcji celu  wnowuym punkcie x1
		} while (X0.y <= X1.y); // je�li prawa stona jest mniejsza to si� zatrzymujmy bo mmay minimum
		if (d > 0) // je�li idziemy w prawo to idziemy 
		{
			p[0] = temp;
			p[1] = m2d(X1.x);
		}
		else
		{
			p[1] = temp;
			p[0] = m2d(X1.x);
		}
		return p; // zwraca przedzia� w kt�rym znaleziono minimum
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution a0(a), b0(b);
		a0.fit_fun(ff, ud1, ud2);
		b0.fit_fun(ff, ud1, ud2);
		solution c(0), d(0); //Punkty wewn�trzne przedzia�u [a , b]
		int k = 1;

		while (GetFib(k) < (b0.x - a0.x) / epsilon) {
			k++;
		}

		c.x = b0.x - GetFib(k - 1) / GetFib(k) * (b0.x - a0.x);
		d.x = a0.x + b0.x - c.x;

		c.fit_fun(ff, ud1, ud2);  // oblczenie wartosci funkcji celu
		d.fit_fun(ff, ud1, ud2);
		Xopt.ud = b - a;
		for (int i = 0; i < k - 3; i++) // punky c i d dziel� przedzia� a i b na 3 cz�sci dlatego k  - 3
		{

			if (c.y < d.y)
			{
				b0 = d; // je�li wartos� funkcji w pounkcie c jest mniejsza ni� w d to minimalna warto�c funkcji znajduje si� w przedziale [a d]
			}
			else
			{
				a0 = c; //a je�li nie to [c b]
			}

			c.x = b0.x - GetFib(k - i - 2) / GetFib(k - i - 1) * (b0.x - a0.x); // w ka�dej itaracji ponownie 
			d.x = a0.x + b0.x - c.x; // obliczamy wsp�rz�dne c.x i d.x na podstaniw zmniejszonego przedzia�u wykorzystuj�c kolejne liczby fibbonaciego
			c.fit_fun(ff, ud1, ud2);
			d.fit_fun(ff, ud1, ud2);
			Xopt.ud.add_row(m2d(b0.x - a0.x));
		}
		Xopt = c;

		Xopt.flag = 0;

		return Xopt;
		
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		double l, m;
		double c = (a + b) / 2.0;
		solution a0(a), b0(b), c0(c), d0(0), di(a);

		a0.fit_fun(ff, ud1, ud2);
		b0.fit_fun(ff, ud1, ud2);
		c0.fit_fun(ff, ud1, ud2);

		// Inicjacja zapisu do pliku
		std::ofstream logFile("lagrange.txt", std::ios::out);
		if (logFile.is_open()) {
			logFile << "Pocz�tek algorytmu:\n";
			logFile << "a0.x = " << a0.x << ", a0.y = " << a0.y << "\n";
			logFile << "b0.x = " << b0.x << ", b0.y = " << b0.y << "\n";
			logFile << "c0.x = " << c0.x << ", c0.y = " << c0.y << "\n\n";
			logFile.close();
		}
		Xopt.ud = b - a;
		do
		{
			Xopt.ud.add_row(m2d(b0.x - a0.x));
			l = m2d(a0.y) * m2d(pow(b0.x, 2) - pow(c0.x, 2)) +
				m2d(b0.y) * m2d(pow(c0.x, 2) - pow(a0.x, 2)) +
				m2d(c0.y) * m2d(pow(a0.x, 2) - pow(b0.x, 2));

			m = m2d(a0.y) * m2d(b0.x - c0.x) +
				m2d(b0.y) * m2d(c0.x - a0.x) +
				m2d(c0.y) * m2d(a0.x - b0.x);

			// Zapis do pliku l i m
			logFile.open("lagrange.txt", std::ios::app);
			if (logFile.is_open()) {
				logFile << "Iteracja " << i << ":\n";
				logFile << "l = " << l << ", m = " << m << "\n";
				logFile.close();
			}

			if (m <= 0)
			{
				Xopt.flag = -1;
				// Zapis b��du do pliku
				logFile.open("lagrange.txt", std::ios::app);
				if (logFile.is_open()) {
					logFile << "B��d: m <= 0.\n\n";
					logFile.close();
				}
				return Xopt;
			}

			di.x = d0.x;
			d0.x = 0.5 * l / m;
			d0.fit_fun(ff, ud1, ud2);
			di.fit_fun(ff, ud1, ud2);

			// Zapis d0 i di do pliku
			logFile.open("lagrange.txt", std::ios::app);
			if (logFile.is_open()) {
				logFile << "d0.x = " << d0.x << ", d0.y = " << d0.y << "\n";
				logFile << "di.x = " << di.x << ", di.y = " << di.y << "\n\n";
				logFile.close();
			}

			if (a0.x < d0.x && d0.x < c0.x)
			{
				if (d0.y < c0.y)
				{
					b0 = c0;
					c0 = d0;
				}
				else
				{
					a0 = d0;
				}
			}
			else
			{
				if (c0.x < d0.x && d0.x < b0.x)
				{
					if (d0.y < c0.y)
					{
						a0 = c0;
						c0 = d0;
					}
					else
					{
						b0 = d0;
					}
				}
				else
				{
					// Zapis b��du do pliku
					logFile.open("lagrange.txt", std::ios::app);
					if (logFile.is_open()) {
						logFile << "B��d: d0.x poza zakresem.\n";
						logFile.close();
					}
					Xopt = d0;
					Xopt.flag = -1;
					return Xopt;
				}
			}

			i = i + 1;
			Xopt.ud.add_row((b0.x - a0.x));

			// Zapis punkt�w do pliku
			logFile.open("lagrange.txt", std::ios::app);
			if (logFile.is_open()) {
				logFile << "Po aktualizacji:\n";
				logFile << "a0.x = " << a0.x << ", a0.y = " << a0.y << "\n";
				logFile << "b0.x = " << b0.x << ", b0.y = " << b0.y << "\n";
				logFile << "c0.x = " << c0.x << ", c0.y = " << c0.y << "\n";
				logFile << "d0.x = " << d0.x << ", d0.y = " << d0.y << "\n\n";
				logFile.close();
			}

			if (solution::f_calls > Nmax)
			{
				// Zapis b��du do pliku
				logFile.open("lagrange.txt", std::ios::app);
				if (logFile.is_open()) {
					logFile << "B��d: Nie znaleziono przedzia�u po " << Nmax << " pr�bach.\n\n";
					logFile.close();
				}
				throw std::runtime_error("Nie znaleziono przedzialu po " + std::to_string(Nmax) + " probach");
			}

		} while ((b0.x - a0.x) >= epsilon && fabs(m2d(d0.x) - m2d(di.x)) >= gamma);

		Xopt = d0;
		Xopt.fit_fun(ff, ud1, ud2);

		// Zapis wyniku do pliku
		logFile.open("lagrange.txt", std::ios::app);
		if (logFile.is_open()) {
			logFile << "Zako�czenie algorytmu:\n";
			logFile << "Xopt.x = " << Xopt.x << ", Xopt.y = " << Xopt.y << "\n\n";
			logFile.close();
		}
		Xopt.flag = 0;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}


solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
