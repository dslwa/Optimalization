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
{//funkcja celu, początkowa warośc zmiennej, krok który określa odlwegłośc między punktem początkowaym a końcowym, alpha rozszerza lub zmienjsza krok z każdym kroku pętli, max liczba wywołań funkcji celu
	try
	{
		double* p = new double[2]{ 0,0 }; //przechowuje nasze przedziały gdzie jest miinmum
		int i = 0;
		double temp = 0; 
		solution X0(x0), X1(x0 + d); //punkt początkowy, obiekty klasy solution x1 = x0 + d
		X0.fit_fun(ff, ud1, ud2); 
		X1.fit_fun(ff, ud1, ud2);// wylicza wartość funkcji celu w tych punktach
		if (X0.y == X1.y)
		{
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x); // Jeśli te same wartosci to zwraca ten przedział bo zakłada, że nie ma w nim znaczacych zmien finkcji celu
			std::cout << solution::f_calls << "";
			solution::clear_calls();
			return p;
		}	
		if (X1.y > X0.y)
		{
			d = -d; // jeśli wartośc funcji w x1 > x0 to funkcja rośnie czyli nie zmierzamy w stronę minimum, zmieniamy kierunek daltego d na - bo chcemy się cofnąć
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2); // liczymy funkcje celu 
			if (X1.y >= X0.y)
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x) - d; // jeśli nadal jest większa to zwracamy ten przedział
				std::cout << solution::f_calls << "";
				solution::clear_calls();
				return p;
			}
		}
		do
		{
			if (solution::f_calls > Nmax) //sprawdza czy liczba wywolan funkcji celu orzekroczyła dozwolną max liczbę
			{
				X0.flag = 0; 
				std::cout << solution::f_calls << "";
				solution::clear_calls();
				return 0;
			}	
			i = i + 1;
			temp = m2d(X0.x);

			X1 = X0.x + pow(alpha, i) * d; //oblicza nowy punkt x1 zwiekszajac krok alfa do poteki i
			X1.fit_fun(ff, ud1, ud2); // liczbymy wartosc funkcji celu  wnowuym punkcie x1
		} while (X0.y <= X1.y); // jeśli prawa stona jest mniejsza to się zatrzymujmy bo mmay minimum
		if (d > 0) // jeśli idziemy w prawo to idziemy 
		{
			p[0] = temp;
			p[1] = m2d(X1.x);
		}
		else
		{
			p[1] = temp;
			p[0] = m2d(X1.x);
		}
		std::cout << solution::f_calls << "";
		solution::clear_calls();
		return p; // zwraca przedział w którym znaleziono minimum
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
		solution Xopt; // zmienna przechowująca wynik
		solution a0(a), b0(b);
		a0.fit_fun(ff, ud1, ud2);
		b0.fit_fun(ff, ud1, ud2);
		solution c(0), d(0); //Punkty wewnętrzne przedziału [a , b]
		int k = 1;

		// Znajdujemy najmniejszą liczbę k spełniającą nierówność φk > (b - a) / ε
		while (GetFib(k) <= (b0.x - a0.x) / epsilon) {
			k++;
		}

		// Upewniamy się, że wartości φk-1 i φk są różne od zera bo nie mozemy dzielić przez 0
		if (GetFib(k - 1) == 0 || GetFib(k) == 0) {
			throw std::runtime_error("Division by zero in Fibonacci calculation");
		}

		// Wyznaczenie punktów wewnętrznych na podstanie długości pfrzedziału liczb fibonaciego
		c.x = b0.x - ((GetFib(k - 1)) / GetFib(k)) * (b0.x - a0.x);
		d.x = a0.x + b0.x - c.x;

		c.fit_fun(ff, ud1, ud2); // oblczenie wartosci funkcji celu
		d.fit_fun(ff, ud1, ud2);


		for (int i = 0; i <= k - 3; i++) // punky c i d dzielą przedział a i b na 3 częsci dlatego k  - 3
		{
		std::cout << "Iteracja " << i << ": c.x = " << m2d(c.x) << ", d.x = " << m2d(d.x) << "\n"; // wypisujemy krokoi poszukiwania minimum
		std::cout << "Wartości funkcji: c.y = " << m2d(c.y) << ", d.y = " << m2d(d.y) << "\n";

			if (c.y < d.y)
			{
				b0 = d; // jeśli wartosć funkcji w pounkcie c jest mniejsza niż w d to minimalna wartośc funkcji znajduje się w przedziale [a d]
			}
			else
			{
				a0 = c; //a jeśli nie to [c b]
			}


			c.x = b0.x - (((GetFib(k - i - 2)) / GetFib(k - i - 1)) * (b0.x - a0.x)); // w każdej itaracji ponownie 
			d.x = a0.x + b0.x - c.x; // obliczamy współrzędne c.x i d.x na podstaniw zmniejszonego przedziału wykorzystując kolejne liczby fibbonaciego

			c.fit_fun(ff, ud1, ud2);
			d.fit_fun(ff, ud1, ud2); // wartości funkcji celu w nowych punktach

			// Przerwanie jeśli różnica pomiędzy a(i) i b(i) jest mniejsza niż epsilon
			if (fabs(m2d(b0.x) - m2d(a0.x)) < epsilon) {
				break;
			}
		}

		// Zwracamy ostateczne rozwiązanie czyli minimum funkcji
		Xopt = c;
		std::cout << solution::f_calls << ",";
		solution::clear_calls();
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


			if (m <= 0)
			{
				Xopt.flag = -1;
				std::cout << solution::f_calls << ",";
				solution::clear_calls();
				return Xopt;
			}

			di.x = d0.x;
			d0.x = 0.5 * l / m;
			d0.fit_fun(ff, ud1, ud2);
			di.fit_fun(ff, ud1, ud2);


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
					Xopt = d0;
					Xopt.flag = -1;
					std::cout << solution::f_calls << ",";
					solution::clear_calls();
					return Xopt;
				}
			}

			i = i + 1;
			Xopt.ud.add_row((b0.x - a0.x));

			if (solution::f_calls > Nmax)
			{
				std::cout << "B³¹d: Nie znaleziono przedzia³u po " << Nmax << " próbach.\n\n";
				throw std::runtime_error("Nie znaleziono przedzialu po " + std::to_string(Nmax) + " probach");
			}

		} while ((b0.x - a0.x) >= epsilon && fabs(m2d(d0.x) - m2d(di.x)) >= gamma);

		Xopt = d0;
		Xopt.fit_fun(ff, ud1, ud2);
		Xopt.flag = 0;
		
		std::cout << solution::f_calls << ",";
		solution::clear_calls();
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
