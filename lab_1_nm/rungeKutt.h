#pragma once
#include<iostream>
#include<vector>
#include<cmath>



class RungeKutt
{
private:
	double testF(double x, double u)
	{
		if (ans == 1)	
		{
			double f = u;
			return f;
		}
		else if (ans == 2)
		{
			double f = (x / (1 + pow(x, 2))) * pow(u, 2) + u - pow(u, 3) * sin(10 * x);
			return f;
		}
	}
	int ans;

protected:
	int size;
	int const MAXSIZE = 100000;
	double hop, e;
	double x;
	double u;
	std::vector<double> rk4(double x, double u, double h)
	{

		double k1 = testF(x, u);
		double k2 = testF(x + 0.5 * h, u + 0.5 * h * k1);
		double k3 = testF(x + 0.5 * h, u + 0.5 * h * k2);
		double k4 = testF(x + 0.5 * h, u + 0.5 * h * k3);

		x = x + h;
		u = u + h * (k1 + 2 * (k2 + k3) + k4) / 6;
		std::vector<double> a = { x, u };
		return a;

	}
	std::vector<double> rk4WithHalfHop(double x, double u, double h)
	{
		return rk4(x, u, h / 2);
	}
	std::vector<double> rk4WithCorrectHop(double x, double u, double h)
	{
		double j;
		j = correctHop(x, u, h);

		return rk4(x, u, j);

	}
	double correctHop(double x, double u, double h)
	{
		auto a = this->rk4(x, u, h);
		auto b = this->rk4WithHalfHop(x, u, h);

		double s = (b[1] - a[1]) / (pow(2, 4) - 1);
		double lec = s * pow(2, 4);

		if (abs(s) > e)
		{
			static int count1 = 0;
			count1++;
			return correctHop(x, u, h / 2);
		}
		else if (abs(s) >= pow(e / 2, 5) && abs(s) <= e)
		{
			return h;
		}
		else if (abs(s) < pow(e / 2, 5))
		{
			static int count2 = 0;
			count2++;
			return 2 * h;
		}
	}
public:

	RungeKutt(int _size = 100, double _h = 0.01, double _e = 0.000001, double x0 = 0, double u0 = 1,int _ans = 1) :size(_size), hop(_h), e(_e), x(x0), u(u0), ans(_ans)
	{
	}
	~RungeKutt() {}
	RungeKutt(const RungeKutt& other)
	{
		this->size = other.size;
		this->x = other.x;
		this->u = other.u;
		this->hop = other.hop;
		this->e = other.e;
		ans = other.ans;
	}
	virtual std::vector<std::vector<double>> rk4()
	{
		std::vector<std::vector<double>> res;
		std::vector<double> vecX, vecU;
		vecX.push_back(this->x);
		vecU.push_back(this->u);

		for (std::size_t i = 0; i < size && i < MAXSIZE; i++)
		{
			auto a = this->rk4(vecX.back(), vecU.back(), hop);
			vecX.push_back(a[0]);
			vecU.push_back(a[1]);
		}
		res.push_back(vecX);
		res.push_back(vecU);
		return res;
	}
	virtual std::vector<std::vector<double>> rk4WithHalfHop()
	{
		std::vector<std::vector<double>> res;
		std::vector<double> vecX, vecU;
		vecX.push_back(this->x);
		vecU.push_back(this->u);

		for (std::size_t i = 0; i < size && i < MAXSIZE; i++)
		{
			auto a = this->rk4WithHalfHop(vecX.back(), vecU.back(), hop);
			vecX.push_back(a[0]);
			vecU.push_back(a[1]);
		}
		res.push_back(vecX);
		res.push_back(vecU);
		return res;
	}
	virtual std::vector<std::vector<double>> rk4WithCorrectHop()
	{
		std::vector<std::vector<double>> res;
		std::vector<double> vecX, vecU;
		vecX.push_back(this->x);
		vecU.push_back(this->u);

		for (std::size_t i = 0; i < size && i < MAXSIZE; i++)
		{
			auto a = this->rk4WithCorrectHop(vecX.back(), vecU.back(), hop);
			vecX.push_back(a[0]);
			vecU.push_back(a[1]);
		}
		res.push_back(vecX);
		res.push_back(vecU);
		return res;
	}
	void printRK4()
	{
		auto res = this->rk4();
		std::cout << "x: ";
		for (int i = 0; i < this->size; i++)
		{
			std::cout<< res[0].at(i)<<" ";
		}
		std::cout << std::endl;
		std::cout << "u: ";
		for (int i = 0; i < this->size; i++)
		{
			
			std::cout << res[1].at(i)<<" ";
			
		}
		std::cout << std::endl;
		std::cout << std::endl;

	}
	void printRKwhs4()
	{
		auto res = this->rk4WithHalfHop();
		std::cout << "x: ";
		for (int i = 0; i < this->size; i++)
		{
			std::cout << res[0].at(i) << " ";
		}
		std::cout << std::endl;
		std::cout << "u: ";
		for (int i = 0; i < this->size; i++)
		{

			std::cout << res[1].at(i) << " ";

		}
		std::cout << std::endl;
		std::cout << std::endl;
	}

	void printRKwcs4()
	{
		auto res = this->rk4WithCorrectHop();
		std::cout << "x: ";
		for (int i = 0; i < this->size; i++)
		{
			std::cout << res[0].at(i) << " ";
		}
		std::cout << std::endl;
		std::cout << "u: ";
		for (int i = 0; i < this->size; i++)
		{

			std::cout << res[1].at(i) << " ";

		}
		std::cout << std::endl;
		std::cout << std::endl;
	}

};

class RungeKuttForSystem:public RungeKutt
{
private:
	double u1;
	double main2F(double x, double u, double u1, int num)
	{
		int a = 3;
		int b = 3;
		if (num == 1)
			return u1;
		else if (num == 2)
			return -a * u1 + b * sin(u);

	}
	std::vector<double>rk4(double x, double u, double u1, double h)
	{
		double k11 = main2F(x, u, u1, 1);
		double k12 = main2F(x, u, u1, 2);

		double k21 = main2F(x + 0.5 * h, u + 0.5 * h * k11, u1 + 0.5 * h * k12, 1);
		double k22 = main2F(x + 0.5 * h, u + 0.5 * h * k11, u1 + 0.5 * h * k12, 2);

		double k31 = main2F(x + 0.5 * h, u + 0.5 * h * k21, u1 + 0.5 * h * k22, 1);
		double k32 = main2F(x + 0.5 * h, u + 0.5 * h * k21, u1 + 0.5 * h * k22, 2);

		double k41 = main2F(x + 0.5 * h, u + 0.5 * h * k31, u1 + 0.5 * h * k32, 1);
		double k42 = main2F(x + 0.5 * h, u + 0.5 * h * k31, u1 + 0.5 * h * k32, 2);


		x = x + h;
		u = u + h * (k11 + 2 * (k21 + k31) + k41) / 6;
		u1 = u1 + h * (k12 + 2 * (k22 + k32) + k42) / 6;
		std::vector<double> a = { x, u,u1 };
		return a;
	}
	std::vector<double> rk4WithHalfHop(double x, double u, double u1, double h)
	{
		return rk4(x, u, u1, h / 2);
	}
	std::vector<double> rk4WithCorrectHop(double x, double u, double u1, double h)
	{
		double j;
		j = correctHop(x, u, u1, h);

		return rk4(x, u, u1, j);
	}
	double correctHop(double x, double u, double u1, double h)
	{
		auto a = this->rk4(x, u, u1, h);
		auto b = this->rk4WithHalfHop(x, u, u1, h);

		double s = (a[2] - b[2]) / (pow(2, 4) - 1);
		double lec = s * pow(2, 4);

		if (abs(s) > e)
		{
			static int count1 = 0;
			count1++;
			return correctHop(x, u, u1, h / 2);
		}
		else if (abs(s) >= pow(e / 2, 5) && abs(s) <= e)
			return h;
		else if (abs(s) < pow(e / 2, 5))
		{
			static int count2 = 0;
			count2++;
			return 2 * h;
		}
	}

public:
	RungeKuttForSystem(int _size , double _h , double _e , double x0, double u0 , double _u1 ):RungeKutt(_size,_h,_e,x0,u0), u1(_u1)
	{}
	~RungeKuttForSystem(){}
	std::vector<std::vector<double>> rk4()
	{
		std::vector<std::vector<double>> res;
		std::vector<double> vecX, vecU, vecU1;
		vecX.push_back(this->x);
		vecU.push_back(this->u);
		vecU1.push_back(this->u1);

		for (size_t i = 0; i < size && i < MAXSIZE; i++)
		{
			auto a = this->rk4(vecX.back(), vecU.back(), vecU1.back(), hop);
			vecX.push_back(a[0]);
			vecU.push_back(a[1]);
			vecU1.push_back(a[2]);
		}
		res.push_back(vecX);
		res.push_back(vecU);
		res.push_back(vecU1);
		return res;
	} 
	std::vector<std::vector<double>> rk4WithHalfHop()
	{
		std::vector<std::vector<double>> res;
		std::vector<double> vecX, vecU, vecU1;
		vecX.push_back(x);
		vecU.push_back(u);
		vecU1.push_back(u1);

		for (size_t i = 0; i < size && i < MAXSIZE; i++)
		{
			auto a = this->rk4WithHalfHop(vecX.back(), vecU.back(), vecU1.back(), hop);
			vecX.push_back(a[0]);
			vecU.push_back(a[1]);
			vecU1.push_back(a[2]);
		}
		res.push_back(vecX);
		res.push_back(vecU);
		res.push_back(vecU1);
		return res;
	}
	std::vector<std::vector<double>> rk4WithCorrectHop()
	{
		std::vector<std::vector<double>> res;
		std::vector<double> vecX, vecU, vecU1;
		vecX.push_back(this->x);
		vecU.push_back(this->u);
		vecU1.push_back(this->u1);

		for (size_t i = 0; i < size && i < MAXSIZE; i++)
		{
			auto a = this->rk4WithCorrectHop(vecX.back(), vecU.back(), vecU1.back(), hop);
			vecX.push_back(a[0]);
			vecU.push_back(a[1]);
			vecU1.push_back(a[2]);
		}
		res.push_back(vecX);
		res.push_back(vecU);
		res.push_back(vecU1);
		return res;
	}

};