//============================================================================
// Name        : Kursusopgave.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

// Headerfiler:

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
using namespace std;

// Konstanter og variable:

#define NMAX 10

// Delprogrammer:
void IndtastEllerHentData(double A[NMAX][NMAX], double y[NMAX], int &n, double K[NMAX][NMAX+1]);
void Del1(double y[NMAX], double A[NMAX][NMAX], double b[NMAX], double c[NMAX], double K[NMAX][NMAX+1], int n, int &bs, double x[NMAX]);
void Del2(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod);
void Del3(double A[NMAX][NMAX], double RodTabel[NMAX], int n, double EigVecs[NMAX][NMAX], double EigVals[NMAX], int &NumbEigSolsFound);

// Hjaelpefunktioner:
	// Generelle hjaelpefunktioner og introduktioner:

void Intro();
void IntroDel1();
void IntroDel2();
void IntroDel3();
void IndhentTotalmatrix(double A[NMAX][NMAX],double y[NMAX],int &n);
void IndtastData(double A[NMAX][NMAX], double y[NMAX], int &n);
void DanTotalMatrix(double K[NMAX][NMAX+1],double A[NMAX][NMAX], double y[NMAX], int n);
void MatrixUdskrivTM(double K[NMAX][NMAX+1], int n);
void MatrixUdskriv(double M[NMAX][NMAX], int n);

	// Del 1:

void Krylovs(double y[NMAX], double A[NMAX][NMAX], double b[NMAX], double c[NMAX], double K[NMAX][NMAX+1], int n);
void Kopicib(double v[NMAX], double b[NMAX], int n);
void IndsaetKrylov(double K[NMAX][NMAX+1], double y[NMAX], int j, int n);
void Produkt(double M[NMAX][NMAX], double v[NMAX], double vprod[NMAX], int n);
void BackwardsSubstitution(double Alb[NMAX][NMAX+1], double x[NMAX], int n);
void Gauss(double TM[NMAX][NMAX+1], int n, int &bs);
void DelvisPivotering(double TM[NMAX][NMAX+1], int j, int n);
void Singularitet();
void UdskrivPolynomium(double koef[NMAX + 1], int n);
void ArrangerKoef(double koef[NMAX], int n);

	// Del 2:

void Horner(int n, double koef[NMAX+1], double x, double &fx, double &dfdx);
void UdskrivHornerTabel(double koef[NMAX], int n);
void IndtastNRData(double &x0, double &eps, int &N);
void NewRapHorner(double x0, double eps, int N, int n, double koef[NMAX+1],double &rod, int &succes);
void HornerDivPol(int &n, double koef[NMAX+1],double rod);

	// Del 3:

void NormerVektor(double v[NMAX],double vnorm[NMAX],double &l,int n);
void VektorPrint(double v[NMAX], int n);
void DanUdvKoefMatInvIt(double A[NMAX][NMAX], double lamst, double yk_1[NMAX], int n, double UdvKoefMat[NMAX][NMAX+1]);
double LgdDifVect(double v1[NMAX], double v2[NMAX], int n);
void IndtastVektor(double v[NMAX], int n);
void IndtastIIData(double &dlamb, double &e1, double &e2, int &N);
void IndsaetInvers(double K[NMAX][NMAX], double y[NMAX], int j, int n);
void InversIteration(double A[NMAX][NMAX], double RodTabel[NMAX], int n, double EigVecs[NMAX][NMAX], double EigVals[NMAX], int &NumbEigSolsFound);

// Main-funktion:

int main() {
	double x0, eps, rod, y[NMAX], A[NMAX][NMAX], b[NMAX],c[NMAX], K[NMAX][NMAX+1], x[NMAX],
	RodTabel[NMAX], EigVecs[NMAX][NMAX], EigVals[NMAX];
	int n, N, succes, bs, NumbEigSolsFound;
	char genkørsel;
	do {
		Intro();

		IntroDel1();
		IndtastEllerHentData(A,y,n,K);

		// Del 1:
		Del1(y,A,b,c,K,n,bs,x);

		// Del 2:
		IntroDel2();
		Del2(N,eps,x0,x,n,succes,RodTabel,rod);

		// Del 3:
		IntroDel3();
		Del3(A,RodTabel,n,EigVecs,EigVals,NumbEigSolsFound);

		cout << endl << endl << "Vil du køre hele programmet forfra? Indtast j/n: "; cin >> genkørsel;
		}
	while (genkørsel == 'j' or genkørsel == 'J');
	return 0;
}

// Delprogrammer:

void IndtastEllerHentData(double A[NMAX][NMAX], double y[NMAX], int &n, double K[NMAX][NMAX+1]){
	int svarData;
	cout << "Vil du: \n1) Indhente data fra fil \n2) Indtaste data manuelt \n";
	do {
		cout << "\nIndtast 1 eller 2: ";
		cin >> svarData;
		}
	while (svarData != 1 && svarData != 2);
	switch (svarData){
		case 1:
			cout << "\nData indhentes fra fil." << endl;
			IndhentTotalmatrix(A,y,n);
			DanTotalMatrix(K,A,y,n);
			MatrixUdskrivTM(K,n);
			break;
		case 2:
			cout << "\nData indtastes manuelt." << endl;
			char retMat = 'j';
			do {
				IndtastData(A,y,n);
				DanTotalMatrix(K,A,y,n);
				MatrixUdskrivTM(K,n);
				cout << "\nVil du rette matricen og/eller vektoren? Indtast j/n: "; cin >> retMat;
			}
			while (retMat == 'j' or retMat == 'J');
			break;
	}
}

void Del1(double y[NMAX], double A[NMAX][NMAX], double b[NMAX], double c[NMAX], double K[NMAX][NMAX+1], int n, int &bs, double x[NMAX]){
	Krylovs(y,A,b,c,K,n);
	Gauss(K,n,bs);
	cout << "\nDen Gauss−eliminerede matrix:\n";
	MatrixUdskrivTM(K,n);
	if(bs==1){
		BackwardsSubstitution(K,x,n);
		cout << "\nDen fundne vektor, der løser det pågældende lineære ligningssystem:\n";
		VektorPrint(x,n);
		cout << "\n\nDet karakteristiske polynomium:";
		ArrangerKoef(x,n);
		UdskrivPolynomium(x,n);
	}
	else{
		Singularitet();
	}
}

void Del2(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod){
	char svarNR, svarTabel;
	do {
		UdskrivHornerTabel(x, n);
		cout << "\nVil du udregne rødder vha. Newton Raphson-metoden? Indtast j/n: "; cin >> svarNR;
		int t = 0;
		if (svarNR=='j' or svarNR=='J'){
			IndtastNRData(x0,eps,N);
			do {
				NewRapHorner(x0,eps,N,n,x,rod,succes);
				RodTabel[t] = rod;
				t = t+1;
				if (succes == 1){
					HornerDivPol(n,x,rod);
					UdskrivPolynomium(x,n);
					if (n==1){
						rod = -(x[1]/x[0]);
						cout << endl << endl << "Den sidste rod er fundet og vises her: rod = " << rod << endl;
						svarNR = 'n';
						RodTabel[t] = rod;
					}
					else {
						svarNR = 'j';
					}
				}
				else {
					cout << "\nNR-iterationen mislykkedes.";
					cout << "\nSe en tabel igen og kontroller, at den viser fortegnsskift";
					cout << "\nEfter konstatering af fortegnsskift: Indtast nye NR-data";
					svarTabel = 'j';
					svarNR = 'n';
				}
			}
			while (svarNR == 'j' or svarNR == 'J');
		}
	}
	while (svarTabel == 'j' or svarTabel == 'J');
}

void Del3(double A[NMAX][NMAX], double RodTabel[NMAX], int n, double EigVecs[NMAX][NMAX], double EigVals[NMAX], int &NumbEigSolsFound){
	InversIteration(A,RodTabel,n,EigVecs,EigVals,NumbEigSolsFound);
	cout << endl << "Egenvektorerne: ";
	MatrixUdskriv(EigVecs,n);
	cout << endl << "Egenvaerdierne: " << endl;
	VektorPrint(EigVals,n);
}


// Hjaelpefunktioner:
	// Generelle hjaelpefunktioner og introduktioner:

void Intro(){
	cout << "Dette program finder egenløsninger for en brugervalgt matrix. "
			"Ved at opstille en \nKrylovmatrix ud fra den oprindelige"
			"matrix og en brugervalgt vektor, findes de normaliserede \nkoefficienter "
			"for det karakteristiske polynomium. \n\nPolynomiet opstilles, og rødderne "
			"for polynomiet, dvs. matricens egenvaerdier, estimeres vha. \nNewton/Raphson/Horner-metoden."
			"\n\nDernaest anvendes invers iteration til at estimere egenvektorerne til de "
			"dertilhørende egenvaerdier. \nTil sidst praesenteres løsningen."
			"\n\nProgrammet er delt i tre underdele , der praesenteres hver for sig i det følgende.\n";
	cout << "\n--------------------------------------------------------------------------------------------------\n";
}

void IntroDel1(){
	cout << "\nNu køres Del 1 af programmet. \n\nHer vil Krylov−matricen opstilles og løses vha. bl.a. "
			"Backwards substitution og Gauss, og de \nnormaliserede koefficienter til det karakteristiske "
			"polynomium vil hermed findes. \n\nDu vil nu blive bedt om at indtaste en vektor og "
			"matrix , eller hente delene fra en fil.\n";
	cout << "\n--------------------------------------------------------------------------------------------------\n\n";
}

void IntroDel2(){
	cout << "\n\n--------------------------------------------------------------------------------------------------\n\n";
	cout << "Programmet går nu videre til Del 2. \n\nHer vil programmet forsøge at udlede rødderne, og "
			"herved egenvaerdierne til det fundne polynomium, \nDette gøres vha. Horner og Newton−Raphson metoderne.\n";
	cout << "\n--------------------------------------------------------------------------------------------------\n";
}

void IntroDel3(){
	cout << "\n--------------------------------------------------------------------------------------------------\n";
	cout << "\nDel 3 af programmet påbegyndes nu. \n\nDenne del har til formål at finde egenvektorerne "
			"med tilhørende estimerede \negenværdier (for at vi kan se, de stemmer overens med de allerede "
			"fundne egenværdier og \napproksimeringen er korrekt). \n\nDette gøres vha. invers iteration .\n";
	cout << "\n--------------------------------------------------------------------------------------------------\n";
}

void IndhentTotalmatrix(double A[NMAX][NMAX],double y[NMAX],int &n){
	ifstream Fil;
	Fil.open("TotalMatrix.txt");
	Fil>>n;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++) Fil>>A[i][j];
		Fil>>y[i];
	}
	Fil.close();
}

void IndtastData(double A[NMAX][NMAX], double b[NMAX], int &n){
	cout << "Indtast antal raekker og søjler, nxn: "; cin >> n;
	cout << "Først indtastes matricen A: \n";
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				cout << "Indtast element i raekke nr. " << i+1 << " og søjle nr. " << j+1 << " her: ";
				cin >> A[i][j];
			}
		}
		cout << "Nu indtastes vektoren b: \n";
		for(int i=0; i<n; i++){
			cout << "Indtast element i raekke nr. " << i+1 << " her: ";
			cin >> b[i];
		}
}

void DanTotalMatrix(double TM[NMAX][NMAX+1],double A[NMAX][NMAX], double b[NMAX], int n){
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			TM[i][j] = A[i][j];
		}
	}
	for(int i=0; i<n; i++){
		TM[i][n] = b[i];
	}
}

void MatrixUdskrivTM(double TM[NMAX][NMAX+1], int n){
	for(int i=0; i<n; i++){
		for(int j=0; j<n+1; j++){
			cout << setw(8) << setprecision(5) << TM[i][j];
		}
	cout << "\n";
	}
}

void MatrixUdskriv(double M[NMAX][NMAX], int n){
	cout << endl;
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			cout << setw(10) << setprecision(4) << M[i][j];
		}
		cout << endl;
	}
}

void VektorPrint(double v[NMAX], int n){
	cout << "[";
	for(int k=0; k<n-1; k++){
		cout << setprecision(14) << v[k] << ",";
	}
	cout << setprecision(14) << v[n-1] << "]";
}

	// Del 1 hjaelpefunktioner:

void Krylovs(double y[NMAX], double A[NMAX][NMAX], double b[NMAX], double c[NMAX], double K[NMAX][NMAX+1], int n){
	int j;
	cout << "\nKrylovmatrix:\n";
	Kopicib(y,b,n);
	IndsaetKrylov(K,b,n-1,n);
	for(j=n-2; j>=0; j--){
		Produkt(A,b,c,n);
		IndsaetKrylov(K,c,j,n);
		Kopicib(c,b,n);
	}
	Produkt(A,b,c,n);
	for(j=0; j<n; j++){
		c[j] = -c[j];
	}
	IndsaetKrylov(K,c,n,n);
	MatrixUdskrivTM(K,n);
}

void Kopicib(double v[NMAX], double b[NMAX], int n){
	for(int i=0; i<n; i++){
		b[i] = v[i];
	}
}

void IndsaetKrylov(double K[NMAX][NMAX+1], double y[NMAX], int j, int n){
	for(int i=0; i<n; i++){
		K[i][j] = y[i];
	}
}

void Produkt(double M[NMAX][NMAX], double v[NMAX], double vprod[NMAX], int n){
	for(int i=0; i<n; i++){
		vprod[i] = 0;
		for(int j=0; j<n; j++){
			vprod[i] = vprod[i]+M[i][j]*v[j];
		}
	}
}

void Gauss(double TM[NMAX][NMAX+1], int n, int &bs){
	double factor, eps=1e-9;
	int j, k;
	bs=1;
	for(j=0;j<=n-2;j++){
		DelvisPivotering(TM,j,n);
		if(fabs(TM[j][j])<eps){
			bs=0;
			break;
		}
		for(int i=j+1; i<=n-1; i++){
			factor=-TM[i][j]/TM[j][j];
			TM[i][j]=0;
			for(k=j+1; k<=n; k++){
				TM[i][k]+=factor*TM[j][k];
			}
		}
	}
	if(fabs(TM[n-1][n-1])<eps){
		bs=0;
	}
}

void BackwardsSubstitution(double Alb[NMAX][NMAX+1], double x[NMAX], int n){
	double sum;
	x[n-1] = Alb[n-1][n]/Alb[n-1][n-1];
	for(int i=n-2; i>=0; i--){
		sum=0;
		for(int j=i+1; j<=n-1; j++){
			sum+=Alb[i][j]*x[j];
		}
		x[i]=(Alb[i][n]-sum)/Alb[i][i];
	}
}

void DelvisPivotering(double TM[NMAX][NMAX+1], int j, int n){
	double max, temp;
	int maxpos, i, k;
	max = fabs(TM[j][j]);
	maxpos = j;
	for(i=j+1; i<=n-1; i++){
		if(fabs(TM[i][j]) > max){
			max = fabs(TM[i][j]);
			maxpos = i;
		}
	}
	for (k=j; k<=n; k++){
		temp = TM[j][k];
		TM[j][k] = TM[maxpos][k];
		TM[maxpos][k] = temp;
	}
}

void Singularitet(){
	cout << "\nUnder Gauss-eliminationen har det vist sig, at A formodes at vaere singulaer (ikke-invertibel). "
			"\nBackwards substitution kan derfor ikke udføres.";
}

void ArrangerKoef(double koef[NMAX], int n){
	double temp;
	for (int k=n; k>=0; k--){
		temp = koef[k-1];
		koef[k] = temp;
	}
	if (n%2 != 0){
		koef[0] = -1;
	}
	else {
		koef[0] = 1;
	}
}

void UdskrivPolynomium(double koef[NMAX+1], int n){
	cout << endl << "P(λ)=";
	if (n>1){
		if (n%2 != 0){
			cout << "-";
		}
		cout << "λ^" << n;
		for (int k=1; k<n-1; k++){
			if (koef[k] > 0){
				cout << "+";
			}
			if (koef[k]!=0){
				cout << koef[k] << "λ^" << n-k;
			}
		}
		if (koef[n-1]!=0){
			if (koef[n-1] > 0) {
				cout << "+";
			}
			cout << koef[n-1] << "λ";
		}
		if (koef[n]!=0){
			if (koef[n] >= 0) {
				cout << "+";
			}
			cout << koef[n];
		}
	}
	else if (n==1){
		if (koef[0]!=0){
			cout << koef[0] << "λ";
		}
		if (koef[n]!=0){
			if (koef[n] > 0){
				cout << "+";
			}
			cout << koef[n];
		}
	}
	else if (n==0){
		cout << koef[n];
	}
}

	// Del 2 hjaelpefunktioner:

void Horner(int n, double koef[NMAX+1], double x, double &fx, double &dfdx){
	int k;
	double b[NMAX+1], c[NMAX];
	b[0]=koef[0];
	for(k=1;k<=n;k++){
		b[k]=koef[k]+x*b[k-1];
	}
	c[0]=b[0];
	for(k=1;k<n;k++){
		c[k]=b[k]+x*c[k-1];
	}
	fx=b[n];
	dfdx=c[n-1];
}

void UdskrivHornerTabel(double koef[NMAX], int n){
	double a, b, x, interval, fx, dfdx;
	int funk_vaerdier=20;
	cout << "\nUdskrivning af tabel med Horners vaerdier:";
	cout << "\nProgrammet udregner " << funk_vaerdier << " funktionsvaerdier af funktionen og den afledte funktion i et givent interval [a;b].";
	cout << "\n\nVaelg startvaerdi for x-intervallet: "; cin >> a;
	cout << "Vaelg slutvaerdi for x-intervallet: "; cin >> b;
	interval = fabs(a-b)/(funk_vaerdier-1);
	x = a;
	cout << endl << setw(30) << "x" << "  |  " << setw(30) << "p(x)" << "  |  " << setw(30) << "p'(x)" << endl
			<< "----------------------------------------------------------------------------------------------------------------" << endl;
	for(int k=0; k<funk_vaerdier; k++){
		Horner(n,koef,x,fx,dfdx);
		cout << setprecision(14) << setw(30) << x << "  |  " << setw(30) << setprecision(14) << fx << "  |  "
				<< setw(30) << setprecision(14) << dfdx << endl;
		x = x+interval;
	}
}

void IndtastNRData(double &x0, double &eps, int &N){
	cout << "\nIndtast startpunkt, x0: "; cin >> x0;
	cout << "Indtast N: "; cin >> N;
	cout << "Indtast epsilon: "; cin >> eps;
}

void NewRapHorner(double x0, double eps, int N, int n, double koef[NMAX+1], double &rod, int &succes){
	double x, fx, dfdx;
	int k;
	k = 0;
	x = x0;
	succes = 0;
	rod = 0;
	do{
		Horner(n,koef,x,fx,dfdx);
		if(dfdx==0.0){
			cout << "\nDa f'(x)=0, kan Newton Raphson ikke bruges og afsluttes derfor.";
			break;
		}
		else{
			x = x-(fx/dfdx);
			k += 1;
		}
	}
	while(fabs(fx)>eps && k<N);
	cout << endl << endl << "Iterationer = " << k;
	if(k<N){
		rod = x;
		succes = 1;
		cout << "\nNR-iteration gav rod = " << setprecision(14) << rod << ", som vil blive elimineret" << endl;
	}
	else{
		cout << "\nBestemmelse af en rod mislykkedes, idet max iterationer er nået.";
	}
}

void HornerDivPol(int &n, double koef[NMAX+1],double rod){
	int k;
	for(k=1; k<=n; k++){
		koef[k] = koef[k]+rod*koef[k-1];
	}
	n = n-1;
}

	// Del 3 hjaelpefunktioner:

void NormerVektor(double v[NMAX],double vnorm[NMAX],double &l,int n){
	double temp[NMAX];
	for (int k=0; k<n; k++){
		temp[k] = v[k];
	}
	for (int i=1; i<n; ++i){
		if (fabs(temp[0]) < fabs(temp[i]))
			temp[0] = temp[i];
	    }
	l = temp[0];
	for(int i=0; i<n; i++){
		vnorm[i] = (1/l)*v[i];
	}
}

void DanUdvKoefMatInvIt(double A[NMAX][NMAX], double lamst, double yk_1[NMAX], int n, double UdvKoefMat[NMAX][NMAX+1]){
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			UdvKoefMat[i][j] = A[i][j];
		}
		UdvKoefMat[i][i] = UdvKoefMat[i][i]-lamst;
		UdvKoefMat[i][n] = yk_1[i];
	}
}

double LgdDifVect(double v1[NMAX], double v2[NMAX], int n){
	double sum=0, vd[NMAX];
	for (int i=0; i<n; i++){
		vd[i] = v1[i]-v2[i];
	}
	for (int i=0; i<n; i++){
		sum += pow(vd[i],2);
	}
	return sqrt(sum);
}

void IndtastVektor(double v[NMAX], int n){
	cout << endl;
	for(int k=0; k<n; k++){
		cout << "Indtast element nr. " << k+1 << " her: ";
		cin >> v[k];
	}
}

void IndtastIIData(double &dlamb, double &eps1, double &eps2, int &N){
	cout << "\nIndtast forskydning, Δλ: "; cin >> dlamb;
	if (dlamb>0){
		dlamb = -fabs(dlamb);
	}
	cout << "Indtast tolerance, E1: "; cin >> eps1;
	cout << "Indtast tolerance, E2: "; cin >> eps2;
	cout << "Indtast max antal iterationer, N: "; cin >> N;
}

void IndsaetInvers(double K[NMAX][NMAX], double y[NMAX], int j, int n){
	for(int i=0; i<n; i++){
		K[i][j] = y[i];
	}
}

void InversIteration(double A[NMAX][NMAX], double RodTabel[NMAX], int n, double EigVecs[NMAX][NMAX], double EigVals[NMAX], int &NumbEigSolsFound){
	double Alb[NMAX][NMAX+1], z[NMAX], yGl[NMAX], yNy[NMAX], LGl, LNy, dlamb, lamst, eps1, eps2, LgdyNyMinyGl;
	int N, i, k, bs;
	char ReIterate;
	IndtastIIData(dlamb,eps1,eps2,N);
	NumbEigSolsFound = 0;
	for(i=1; i<=n; i++){
		do {
			IndtastVektor(z,n);
			k = 0;
			NormerVektor(z,yNy,LNy,n);
			lamst = RodTabel[i-1] + dlamb;
			do {
				k = k+1;
				LGl = LNy;
				Kopicib(yNy,yGl,n);
				DanUdvKoefMatInvIt(A,lamst,yNy,n,Alb);
				Gauss(Alb,n,bs);
				if (bs == 0){
					Singularitet();
					cout << "\nBemaerk: 0-tolerancen i Gauss-funktionens singularitetstest har indflydelse her!";
					k = N;
				}
				else {
					BackwardsSubstitution(Alb,z,n);
					NormerVektor(z,yNy,LNy,n);
				}
				LgdyNyMinyGl = LgdDifVect(yNy,yGl,n);
			}
			while (fabs(LNy-LGl)>eps1 && LgdyNyMinyGl>eps2 && k<N);
			if (k<N){
				NumbEigSolsFound = NumbEigSolsFound + 1;
				for (int j = 0; j < n; j++) {
					if (fabs(yNy[j]) < eps1 or fabs(yNy[j]) < eps2) {
						yNy[j] = 0.0;
					}
				}
				IndsaetInvers(EigVecs,yNy,NumbEigSolsFound-1,n);
				EigVals[NumbEigSolsFound-1] = lamst + 1/fabs(LNy);
				cout << endl << "Runde nr. " << i << ": " << endl;
				cout << "Egenvaerdi: " << setprecision(14) << EigVals[NumbEigSolsFound-1] << endl;
				cout << "Egenvektor: " << endl;
				VektorPrint(yNy,n);
				cout << endl << "Antal gennemløb: " << k << endl;
				cout << endl << "---------------------------------------------------------\n";
				ReIterate = 'n';
			}
			else {
				cout << "\nInvers Iteration ud fra rod nr." << i << " mislykkedes";
				cout << "\nTast j hvis du ønsker at vaelge nye parametre og udføre en ny iteration.";
				cout << "\nTast n hvis du ønsker at droppe bestemmelsen af egenvektor nr." << i;
				cout << "\nTast j eller n her: "; cin >> ReIterate;
				if (ReIterate == 'j' or ReIterate == 'J') {
					IndtastIIData(dlamb,eps1,eps2,N);
				}
			}
		}
		while (ReIterate == 'j' or ReIterate == 'J');
	}
}
