#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <functional>
#include <string>

#define mp make_pair

#define PI acos(-1.0)


int N = 100;
using namespace std;

const double EPS = 1e-9;
const double INF = 1e9;

int func_id = 0;

struct Individ {
	int kolv_genes;
	vector <double> genes;

	Individ(int _kolv_genes, double mn = 0.0, double mx = 1.0) {
		kolv_genes = _kolv_genes;
		genes.resize(kolv_genes);
		for (int i = 0; i < kolv_genes; ++i) {
			genes[i] = mn + ((double)rand() / RAND_MAX) * (mx - mn);
		}
	}

	Individ(const vector <double>& _genes) {
		genes = _genes;
		kolv_genes = genes.size();
	}

	Individ() {}

};

double dist(const Individ& indv1, const Individ& indv2) {
	double sum = 0.0;
	for (int i = 0; i < indv1.genes.size(); ++i) {
		double dlt = indv1.genes[i] - indv2.genes[i];
		sum += dlt * dlt;
	}
	return sqrt(sum);
}

double getDelta(int m) {
	double sum = 0.0;
	for (int i = 0; i < m; ++i) {
		double ai = (rand() % m == 0 ? 1.0 : 0.0);
		sum += ai / (double)(1ll << (i + 1));
	}
	return sum;
}

void Mutate(Individ indv, double mn = 0.0, double mx = 1.0, double mmin = 0, double mmax = 1, int m = 10) {
	double alpha = 0.5 * (mx - mn);
	for (int i_gen = 0; i_gen < indv.kolv_genes; ++i_gen) {
		if (rand() % 2) {
			alpha *= -1.0;
		}
		indv.genes[i_gen] += alpha * getDelta(m);
		if (indv.genes[i_gen] < mmin)
			indv.genes[i_gen] = mmin;
		if (indv.genes[i_gen] > mmax)
			indv.genes[i_gen] = mmax;
	}
}

Individ Recombine(const Individ& indv1, const Individ& indv2, double mn = 0.0, double mx = 1.0, double d = 0.25) {
	d = fabs(d);

	Individ child = Individ(indv1.kolv_genes);
	double alpha;
	for (int i_gen = 0; i_gen < indv1.kolv_genes; ++i_gen) {
		alpha = -d + ((double)rand() / RAND_MAX) * (2.0 * d + 1.0);
		child.genes[i_gen] = indv1.genes[i_gen] + alpha * (indv2.genes[i_gen] - indv1.genes[i_gen]);

		if (child.genes[i_gen] < mn)
			child.genes[i_gen] = mn;
		if (child.genes[i_gen] > mx)
			child.genes[i_gen] = mx;

	}
	return child;
}

double Rastrigin_function(vector <double>& x) { // min = 0 ; (0, 0, ... 0)
	double sum = 0.0;
	for (int i = 0; i < x.size(); ++i) {
		sum += x[i] * x[i] - 10 * cos(2 * PI*x[i]);
	}
	return 10 * x.size() + sum;
}

double McCormick_function(vector <double>& x) { // min = 0 ; (0, 0, ... 0)
	double s1 = sin(x[0] + x[1]);
	double s2 = (x[0] - x[1]);
	return s1 + s2*s2 - 1.5*x[0] + 2.5*x[1] + 1;
}

double Booth_function(vector <double>& x) { // min = 0 ; (1, 3)
	double sum1 = x[0] + 2.0 * x[1] - 7.0;
	double sum2 = 2.0 * x[0] + x[1] - 5.0;
	return sum1 * sum1 + sum2 * sum2;
}

double Func3(vector <double>& x) { // min = 0 ; (1, 1)
	double s1 = sin(3.0 * PI * x[0]);
	double s2 = sin(3.0 * PI * x[1]);
	double s3 = sin(2.0 * PI * x[1]);
	return s1 * s1 + (x[0] - 1.0) * (x[0] - 1.0) * (1.0 + s2 * s2) +
		(x[1] - 1.0) * (x[1] - 1.0) * (1.0 + s3 * s3);
}

double Func0(vector <double>& x) {
	double sum = 0.0;
	for (int i = 0; i < x.size(); ++i) {
		sum += x[i] * x[i];
	}
	return sum;
}

double Func4(vector<double>& x)
{
	double s1 = (double)sqrt(abs(x[1] - 0.01 * pow(x[0], 2)));
	double s2 = (double)abs(x[0] + 10);
	return 100 * s1 + 0.01*s2;
		//-cos(x[0])*cos(x[1])*exp(-((x[0] - PI)*(x[0] - PI) + (x[1] - PI)*(x[1] - PI)));
}

double Func5(vector<double>& x)
{
	return 2 * x[0] * x[0] - 1.05*x[0] * x[0] * x[0] * x[0] + (x[0] * x[0] * x[0] * x[0] * x[0] * x[0]) / 6 + x[0] * x[1] + x[1] * x[1];
}

double Func6(vector<double>& x)
{
	return  -(x[1] + 47) * sin(sqrt(abs((double)(x[0] / 2) + (x[1] + 47)))) - x[0] * sin(sqrt(abs(x[0] - (x[1] + 47))));
}

double func(vector <double>& x) {
	if (func_id == 0) {
		return Func4(x);
	}
	else if (func_id == 1) {
		return Func5(x);
	}
	else if (func_id == 2) {
		return Func6(x);
	}
	else if (func_id == 3) {
		return Func0(x);
	}

}
double derivFunc(vector <double>& x, int ind, double delta = 1e-6) {
	if (fabs(delta) < EPS) {
		delta = 1e-6;
	}
	vector <double> d_x = x;
	d_x[ind] += delta;
	return (func(d_x) - func(x)) / delta;
}

void Optimize(Individ& indv, double mmin, double mmax, double alpha = 1e-3) {
	double f_delta = INF;
	int cnt = 0;
	int MAX_CNT = 300;
	while (f_delta > EPS && cnt < MAX_CNT) {
		cnt++;
		double prev_f = func(indv.genes);
		for (int i = 0; i < indv.genes.size(); ++i) {
			double deriv = derivFunc(indv.genes, i);
			indv.genes[i] -= alpha * deriv;
			if (indv.genes[i] < mmin)
				indv.genes[i] = mmin;
			if (indv.genes[i] > mmax)
				indv.genes[i] = mmax;
		}
		double cur_f = func(indv.genes);
		f_delta = fabs(cur_f - prev_f);
	}
}

void getOrder(vector <Individ>& population, vector <int>& result_inds) {
	vector <pair <double, int> > ord(population.size());
	for (int i = 0; i < population.size(); ++i) {
		ord[i] = mp(func(population[i].genes), i);
	}
	sort(ord.begin(), ord.end());
	result_inds.resize(ord.size());
	for (int i = 0; i < ord.size(); ++i) {
		result_inds[i] = ord[i].second;
	}
}

void Selection(vector <Individ>& population, int N, double p = 0.4) {
	vector <pair <double, int> > ord(population.size());
	double sum_f = 0.0;
	double mx_val = 1e6;
	for (int i = 0; i < population.size(); ++i) {
		ord[i] = mp(func(population[i].genes), i);
		sum_f += mx_val - ord[i].first;
	}
	vector <double> borders(population.size());
	borders[0] = (mx_val - ord[0].first) / sum_f;
	for (int i = 1; i < population.size(); ++i) {
		borders[i] = borders[i - 1] + (mx_val - ord[i].first) / sum_f;
	}
	vector <Individ> result;
	result.reserve(N);
	for (int i = 0; i < N; ++i) {
		double rnd = (double)rand() / RAND_MAX;
		int ind = lower_bound(borders.begin(), borders.end(), rnd) - borders.begin();
		result.push_back(population[ind]);
	}
	population = result;

}

double MaxDist(const vector <Individ>& population) {
	double mx_dist = -INF;
	for (int i = 0; i < population.size(); ++i) {
		for (int j = i; j < population.size(); ++j) {
			mx_dist = max(mx_dist, dist(population[i], population[j]));
		}
	}
	return mx_dist;
}

int main() {

	setlocale(LC_ALL, "Russian");
	srand(time(NULL));

	
	int DIM;
	cin >> DIM;
	double mn, mx;
	cin >> mn >> mx;
	cin >> func_id;

	int cnt = 0;
	vector <Individ> population(N);
	for (int i = 0; i < N; ++i) {
		population[i] = Individ(DIM, mn, mx);
	}
	
	double mx_dist = INF;
	int MAX_EPOCH = 5000;
	while (mx_dist > EPS && cnt < MAX_EPOCH) {
		cnt++;
		for (int i = 0; i < population.size(); ++i) {
			Optimize(population[i], mn, mx);
		}
		vector <int> ord_inds;
		getOrder(population, ord_inds);
		int kolv = ord_inds.size() / 2;
		for (int i = 0; i < kolv; ++i) {
			int j = rand() % kolv;
			population.push_back(Recombine(population[ord_inds[i]], population[ord_inds[j]], mn, mx));
		}
		for (int i = 0; i < population.size(); ++i) {
			Mutate(population[i], 0.001, 0.1, mn, mx);
		}

		vector <Individ> new_population;
		new_population.reserve(N);
		getOrder(population, ord_inds);
		for (int i = 0; i < N; ++i) {
			new_population.push_back(population[ord_inds[i]]);
		}
		population = new_population;
		mx_dist = MaxDist(population);
		population[0].genes[0] = -10;
		population[0].genes[1] = 1;
	}
	cout << cnt << endl;
	
	vector <pair <double, int> > ord(population.size());
	for (int i = 0; i < population.size(); ++i) {
		ord[i] = mp(func(population[i].genes), i);
	}
	sort(ord.begin(), ord.end());
	int ind = ord[0].second;

	cout << fixed << setprecision(6);
	cout << func(population[ind].genes) << endl;
	for (int i = 0; i < population[ind].genes.size(); ++i) {
		cout << population[ind].genes[i] << ' ';
	}
	cout << endl;

	system("pause");

	return 0;
}