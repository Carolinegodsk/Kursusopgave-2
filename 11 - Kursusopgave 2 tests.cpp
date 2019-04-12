//============================================================================
// Name        : 11.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
using namespace std;

#define NMAX 10


// Delprogrammer:
void IndtastEllerHentData(double A[NMAX][NMAX], double y[NMAX], int &n, double K[NMAX][NMAX+1]);
void Del1(double y[NMAX], double A[NMAX][NMAX], double b[NMAX], double c[NMAX], double K[NMAX][NMAX+1], int n, int &bs, double x[NMAX]);
void Del2(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod);
void Del3(double A[NMAX][NMAX], double RodTabel[NMAX], int n, double EigVecs[NMAX][NMAX], double EigVals[NMAX], int &NumbEigSolsFound);

// Hjælpefunktioner:

void IndhentTotalmatrix(double A[NMAX][NMAX],double y[NMAX],int &n);
void IndtastData(double A[NMAX][NMAX], double y[NMAX], int &n);
void DanTotalMatrix(double K[NMAX][NMAX+1],double A[NMAX][NMAX], double y[NMAX], int n);
void MatrixUdskrivTM(double K[NMAX][NMAX+1], int n);
void MatrixUdskriv(double M[NMAX][NMAX], int n);

void Krylows(double y[NMAX], double A[NMAX][NMAX], double b[NMAX], double c[NMAX], double K[NMAX][NMAX+1], int n);
void Kopicib(double v[NMAX], double b[NMAX], int n);
void IndsaetKrylow(double K[NMAX][NMAX+1], double y[NMAX], int j, int n);
void Produkt(double M[NMAX][NMAX], double v[NMAX], double vprod[NMAX], int n);

void BackwardsSubstitution(double Alb[NMAX][NMAX+1], double x[NMAX], int n);
void Gauss(double TM[NMAX][NMAX+1], int n, int &bs);
void DelvisPivotering(double TM[NMAX][NMAX+1], int j, int n);
void Singularitet();

void UdskrivPolynomium(double koef[NMAX + 1], int n);
void ArrangerKoef(double koef[NMAX], int n);

void Horner(int n, double koef[NMAX+1], double x, double &fx, double &dfdx);
void UdskrivHornerTabel(double koef[NMAX], int n);
void IndtastNRData(double &x0, double &eps, int &N);
void NewRapHorner(double x0, double eps, int N, int n, double koef[NMAX+1],double &rod, int &succes);
void HornerDivPol(int &n, double koef[NMAX+1],double rod);

void NormerVektor(double v[NMAX],double vnorm[NMAX],double &l,int n);
void VektorPrint(double v[NMAX], int n);
void DanUdvKoefMatInvIt(double A[NMAX][NMAX], double lamst, double yk_1[NMAX], int n, double UdvKoefMat[NMAX][NMAX+1]);
double LgdDifVect(double v1[NMAX], double v2[NMAX], int n);
void IndtastVektor(double v[NMAX], int n);

void IndtastIIData(double &dlamb, double &e1, double &e2, int &N);
void IndsaetInvers(double K[NMAX][NMAX], double y[NMAX], int j, int n);
void InversIteration(double A[NMAX][NMAX], double RodTabel[NMAX], int n, double EigVecs[NMAX][NMAX], double EigVals[NMAX], int &NumbEigSolsFound);


// Testfunktioner

void IndhentTestMatrix(double A[NMAX][NMAX],int &n);
void IndhentTestEigVals(double EigVals[NMAX],int &n);
void DanEigValsMat(double D[NMAX][NMAX], int n, double EigVals[NMAX]);
void Transponering(double M[NMAX][NMAX], double MT[NMAX][NMAX], int n);
void MatrixProdukt(double M1[NMAX][NMAX], double M2[NMAX][NMAX], double M3[NMAX][NMAX], int n, int m, int q);
void InverseDiagonalMat(double M3[NMAX][NMAX], int n);
void TestMatrice(double A[NMAX][NMAX], int n);

void IndhentRegulærTestTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n);
void IndhentSingulærTestTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n);
void IndhentNærsingulærTestTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n);

void TestRegulæreSystemer(double A[NMAX][NMAX], double y[NMAX], double TM[NMAX][NMAX+1], int n, int &bs);
void TestSingulæreSystemer(double A[NMAX][NMAX], double y[NMAX], double TM[NMAX][NMAX+1], int n, int &bs);
void TestNærsingulæreSystemer(double A[NMAX][NMAX], double y[NMAX], double TM[NMAX][NMAX+1], int n, int &bs);

void TestToleranceNRHorner(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod);
void TestStartpunktNRHorner(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod);
void TestPræcisionNRHorner(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod);



int main() {
	double koef[NMAX+1], x0, eps, rod, y[NMAX], A[NMAX][NMAX], TM[NMAX][NMAX+1], b[NMAX],
		c[NMAX], K[NMAX][NMAX+1], x[NMAX], v[NMAX], vnorm[NMAX], l, lamst, yk_1[NMAX], UdvKoefMat[NMAX][NMAX+1],
		v1[NMAX], v2[NMAX], vd[NMAX], lgd, dlamb, e1, e2, RodTabel[NMAX], EigVecs[NMAX][NMAX], EigVals[NMAX], D[NMAX][NMAX], AT[NMAX][NMAX],
		P[NMAX][NMAX], Q1[NMAX][NMAX], Q2[NMAX][NMAX], Q3[NMAX][NMAX], M[NMAX][NMAX];
		int n, N, succes, svarData, j, bs, NumbEigSolsFound;
		char genkørsel, svarTabel, svarNR;
//	do {
//		// Lav intro
//
//		// Start
//		IndtastEllerHentData(A,y,n,K);
//
//		// Del 1 - Krylows, Gauss, BS, Udskrivning af Det Karakteristiske Polynomium
//		Del1(y,A,b,c,K,n,bs,x);
//
//		// Del 2 - Horner og NR
//		Del2(N,eps,x0,x,n,succes,RodTabel,rod);
//
//		// Del 3 - Invers iteration:
//		Del3(A,RodTabel,n,EigVecs,EigVals,NumbEigSolsFound);
//
//		cout << "\nVil du køre hele programmet forfra? Indtast j/n: "; cin >> genkørsel;
//		}
//	while (genkørsel == 'j' or genkørsel == 'J');

	// Opgave 1.2:
//	TestMatrice(A,n);
//
//	// Opgave 2.1:
//	TestRegulæreSystemer(A,y,TM,n,bs);
//
//	// Opgave 2.2:
//	TestSingulæreSystemer(A,y,TM,n,bs);
//
//	// Opgave 2.3:
//	TestNærsingulæreSystemer(A,y,TM,n,bs);

//	TestToleranceNRHorner(N,eps,x0,x,n,succes,RodTabel,rod);
//	TestStartpunktNRHorner(N,eps,x0,x,n,succes,RodTabel,rod);
//	TestPræcisionNRHorner(N,eps,x0,x,n,succes,RodTabel,rod);

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
			cout << "Data indhentes fra fil.\n";
			IndhentTotalmatrix(A,y,n);
			DanTotalMatrix(K,A,y,n);
			MatrixUdskrivTM(K,n);
			break;
		case 2:
			cout << "Data indtastes manuelt.\n";
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
	Krylows(y,A,b,c,K,n);
	Gauss(K,n,bs);
	if(bs==1){
		BackwardsSubstitution(K,x,n);
		cout << "\nDet karakteristiske polynomium:";
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
		int t = 0; //tæller til RodTabel-array
		if (svarNR=='j' or svarNR=='J'){
			IndtastNRData(x0,eps,N);
			do {
				NewRapHorner(x0,eps,N,n,x,rod,succes);
				RodTabel[t] = rod; 	//Rodværdien gemmes i array RET EVT. TIL
				t = t+1;
				if (succes == 1){
					HornerDivPol(n,x,rod);
					UdskrivPolynomium(x,n);
					if (n==1){
						rod = -(x[1]/x[0]);
						cout << endl << endl << "Den sidste rod er fundet og vises her: rod = " << rod << endl;
						svarNR = 'n';
						RodTabel[t] = rod; 	//Rodværdien gemmes i array RET EVT. TIL
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
	cout << endl << "Egenværdierne: " << endl;
	VektorPrint(EigVals,n);
}


// Hjælpefunktioner:

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

void Krylows(double y[NMAX], double A[NMAX][NMAX], double b[NMAX], double c[NMAX], double K[NMAX][NMAX+1], int n){
	int j;
	cout << "\nBeregning vha. Krylows:";
	Kopicib(y,b,n);
	IndsaetKrylow(K,b,n-1,n);
	for(j=n-2; j>=0; j--){
		Produkt(A,b,c,n);
		IndsaetKrylow(K,c,j,n);
		Kopicib(c,b,n);
	}
	Produkt(A,b,c,n);
	for(j=0; j<n; j++){
		c[j] = -c[j];
	}
	IndsaetKrylow(K,c,n,n);
	MatrixUdskrivTM(K,n);
}

void Kopicib(double v[NMAX], double b[NMAX], int n){
	for(int i=0; i<n; i++){
		b[i] = v[i];
	}
}

void IndsaetKrylow(double K[NMAX][NMAX+1], double y[NMAX], int j, int n){
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
	cout << "\nUnder Gauss-eliminationen har det vist sig, at A formodes at vaere singulaer. \nBackwards substitution kan derfor ikke udføres.";
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

void UdskrivPolynomium(double koef[NMAX+1], int n) {
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
	int funk_værdier=20;
	cout << "\n\nUdskrivning af tabel med Horners værdier:";
	cout << "\nProgrammet udregner " << funk_værdier << " funktionsværdier af funktionen og den afledte funktion i et givent interval [a;b].";
	cout << "\n\nVælg startværdi for x-intervallet: "; cin >> a;
	cout << "Vælg slutværdi for x-intervallet: "; cin >> b;
	interval = fabs(a-b)/(funk_værdier-1);
	x = a;
	cout << endl << setw(8) << "x" << "  |  " << setw(20) << "p(x)" << "  |  " << setw(20) << "p'(x)" << endl << "------------------------------------------------------------------------------------------" << endl;
	for(int k=0; k<funk_værdier; k++){
		Horner(n,koef,x,fx,dfdx);
		cout << setprecision(4) << setw(8) << x << "  |  " << setw(20) << setprecision(4) << fx << "  |  " << setw(20) << setprecision(4) << dfdx << endl;
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
	int k, j;
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
	cout << endl << endl << "Iterationer = " << k;		//skal slettes
	if(k<N){
		rod = x;
		succes = 1;
		cout << "\nNR-iteration gav rod = " << rod << ", som vil blive elimineret" << endl;
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

void VektorPrint(double v[NMAX], int n){
	cout << "[";
	for(int k=0; k<n-1; k++){
		cout << v[k] << ",";
	}
	cout << v[n-1] << "]";
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
	for(int k=0; k<n; k++){
		cout << "Indtast element nr. " << k+1 << " her: ";
		cin >> v[k];
	}
}

void IndtastIIData(double &dlamb, double &eps1, double &eps2, int &N){
	cout << "\nIndtast forskydning, Δλ: "; cin >> dlamb;
	//kontrollér input, sikr at dlamp er negativ:
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
				cout << "Egenværdi: " << EigVals[NumbEigSolsFound-1] << endl;
				cout << "Egenvektor: " << endl;
				VektorPrint(yNy,n);
				cout << endl << "Antal gennemløb: " << k << endl << endl;
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


// Testfunktioner

//Opgave 1.2
void TestMatrice(double A[NMAX][NMAX], int n){
	double AT[NMAX][NMAX], D[NMAX][NMAX], Q1[NMAX][NMAX], Q2[NMAX][NMAX], Q3[NMAX][NMAX], M[NMAX][NMAX], EigVals[NMAX];
	IndhentTestMatrix(A,n);
	cout << "Indhentet A matrix: ";
	MatrixUdskriv(A,n);
	IndhentTestEigVals(EigVals,n); 				// Her indtastes/indlæses de ønskede egenværdierne og A.
	cout << endl << "Indhentede egenværdier: " << endl;
	VektorPrint(EigVals,n);
	Transponering(A,AT,n); 						// Transponerer matricen A.
	cout << endl << endl << "A^T matrix: ";
	MatrixUdskriv(AT,n);
	MatrixProdukt(AT,A,Q1,n,n,n);				// Tager matrixproduktet af AT & A.
	DanEigValsMat(D,n,EigVals); 				// Danner Matricen med egenværdier i diagonalen.
	cout << endl << "D matrix: ";
	MatrixUdskriv(D,n);
	InverseDiagonalMat(Q1,n); 					// Danner Q1 = (AT *A)-1
	cout << endl << "Invers matrix Q1: ";
	MatrixUdskriv(Q1,n);
	MatrixProdukt(A,Q1,Q2,n,n,n); 				// Danner Q2 = A*Q1
	cout << endl << "Matrix Q2: ";
	MatrixUdskriv(Q2,n);
	MatrixProdukt(D,AT,Q3,n,n,n); 				// Danner Q3 = D*AT
	cout << endl << "Matrix Q3: ";
	MatrixUdskriv(Q3,n);
	MatrixProdukt(Q2,Q3,M,n,n,n); 				// Danner M = Q2*Q3
	cout << endl << "Matrix M: ";
	MatrixUdskriv(M,n);
}

void IndhentTestMatrix(double A[NMAX][NMAX],int &n){
	ifstream Fil;
	Fil.open("Matrix.txt");
	Fil>>n;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++) Fil>>A[i][j];
	}
	Fil.close();
}

void IndhentTestEigVals(double EigVals[NMAX],int &n){
	ifstream Fil;
	Fil.open("EigVals.txt");
	Fil>>n;
	for (int i=0;i<n;i++){
		Fil>>EigVals[i];
	}
	Fil.close();
}

void DanEigValsMat(double D[NMAX][NMAX], int n, double EigVals[NMAX]){
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			D[i][j] = 0;
		}
	}
	for(int i=0; i<n; i++){
		D[i][i] = EigVals[i];
	}
}

void Transponering(double M[NMAX][NMAX], double MT[NMAX][NMAX], int n){
	for (int i=0; i<=n-1; i++){
		for (int j=0; j<=n-1; j++){
			MT[j][i] = M[i][j];
		}
	}
}

void MatrixProdukt(double M1[NMAX][NMAX], double M2[NMAX][NMAX], double M3[NMAX][NMAX], int n, int m, int q){
	int i, j, k;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			M3[i][j] = 0; 			// Initialisering
			for(k=0; k<q; k++) {
				M3[i][j] += M1[i][k]*M2[k][j];
			}
		}
	}
}

void InverseDiagonalMat(double M3[NMAX][NMAX], int n){
	for(int i=0; i<n; i++){
		M3[i][i] = 1/M3[i][i];
		if (fabs(M3[i][i])<0.001){
			M3[i][i] = 0.0;
		}
	}
}

//Opgave 2.1
void IndhentRegulærTestTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n){
	ifstream Fil;
	Fil.open("Test_RegulærTotalMatrix.txt");
	Fil>>n;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++) Fil>>A[i][j];
		Fil>>y[i];
	}
	Fil.close();
}

void TestRegulæreSystemer(double A[NMAX][NMAX], double y[NMAX], double TM[NMAX][NMAX+1], int n, int &bs){
	double x[NMAX];
	cout << "Test af regulære matricer:";
	cout << endl << "Her er der fundet en vilkårlig regulær matrix til testen.";
	IndhentRegulærTestTotalmatrix(A,y,n);
	DanTotalMatrix(TM,A,y,n);
	cout << endl << "Indhentet matrix: " << endl;
	MatrixUdskrivTM(TM,n);
	Gauss(TM,n,bs);
	if(bs==1){
		BackwardsSubstitution(TM,x,n);
		for(int i=0; i<n; i++){
			cout << endl << "x: " << x[i];
		}
		cout << "\nDet karakteristiske polynomium:";
		ArrangerKoef(x,n);
		UdskrivPolynomium(x,n);
	}
	else{
		Singularitet();
	}
}

//Opgave 2.2
void IndhentSingulærTestTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n){
	ifstream Fil;
	Fil.open("Test_SingulærTotalMatrix.txt");
	Fil>>n;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++) Fil>>A[i][j];
		Fil>>y[i];
	}
	Fil.close();
}

void TestSingulæreSystemer(double A[NMAX][NMAX], double y[NMAX], double TM[NMAX][NMAX+1], int n, int &bs){
	double x[NMAX];
	cout << endl << endl << "Test af singulære matricer:";
	cout << endl << "Her er der fundet en vilkårlig singulær matrix til testen.";
	IndhentSingulærTestTotalmatrix(A,y,n);
	DanTotalMatrix(TM,A,y,n);
	cout << endl << "Indhentet matrix: " << endl;
	MatrixUdskrivTM(TM,n);
	Gauss(TM,n,bs);
	if(bs==1){
		BackwardsSubstitution(TM,x,n);
		cout << "\nDet karakteristiske polynomium:";
		ArrangerKoef(x,n);
		UdskrivPolynomium(x,n);
	}
	else{
		Singularitet();
	}
}

//Opgave 2.3

void IndhentNærsingulærTestTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n){
	ifstream Fil;
	Fil.open("Test_NærsingulærTotalMatrix.txt");
	Fil>>n;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++) Fil>>A[i][j];
		Fil>>y[i];
	}
	Fil.close();
}

void TestNærsingulæreSystemer(double A[NMAX][NMAX], double y[NMAX], double TM[NMAX][NMAX+1], int n, int &bs){
	double x[NMAX];
	cout << endl << endl << "Test af nærsingulære matricer:";
	cout << endl << "Her er der fundet en vilkårlig nærsingulær matrix til testen.";
	IndhentNærsingulærTestTotalmatrix(A,y,n);
	DanTotalMatrix(TM,A,y,n);
	cout << endl << "Indhentet matrix: " << endl;
	MatrixUdskrivTM(TM,n);
	Gauss(TM,n,bs);
	if(bs==1){
		BackwardsSubstitution(TM,x,n);
		cout << "\nDet karakteristiske polynomium:";
		ArrangerKoef(x,n);
		UdskrivPolynomium(x,n);
	}
	else{
		Singularitet();
	}
}

//Opgave 3.1
void TestToleranceNRHorner(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod){
	char svarNR;
	eps = 1;
	for(int i=0; i<5; i++){
		n = 4;
		x[0] = -10;
		x[1] = 35;
		x[2] = -50;
		x[3] = 24;
		ArrangerKoef(x,n);

		cout << endl << "Runde nr. " << i+1 << endl;
		cout << endl << "Epsilon er nu: " << eps;

		cout << "\nUdregner rødder vha. Newton Raphson-metoden: ";
		svarNR = 'j';
		int t = 0;
		if (svarNR=='j' or svarNR=='J'){
			x0 = 0.5;
			cout << endl << endl << "Startpunkt, x0: " << x0;
			N = 1000;
			cout << endl << "Max antal iterationer, N: " << N;
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
					svarNR = 'n';
				}
			}
			while (svarNR == 'j' or svarNR == 'J');
		}
		eps = eps/10;
		cout << endl << "-----------------------------------------------------" << endl;
	}
}

//Opgave 3.2
void TestStartpunktNRHorner(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod){
	char svarNR;
	x0 = 0.1;
	for(int i=0; i<5; i++){
		n = 4;
		x[0] = -10;
		x[1] = 35;
		x[2] = -50;
		x[3] = 24;
		ArrangerKoef(x,n);

		cout << endl << "Runde nr. " << i+1 << endl;
		cout << endl << "x0 er nu: " << x0;

		cout << "\nUdregner rødder vha. Newton Raphson-metoden: ";
		svarNR = 'j';
		int t = 0;
		if (svarNR=='j' or svarNR=='J'){
			eps = 0.001;
			cout << endl << endl << "Epsilon, eps: " << eps;
			N = 1000;
			cout << endl << "Max antal iterationer, N: " << N;
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
					svarNR = 'n';
				}
			}
			while (svarNR == 'j' or svarNR == 'J');
		}
		x0 = x0+1;
		cout << endl << "-----------------------------------------------------" << endl;
	}
}

//Opgave 3.3
void TestPræcisionNRHorner(int &N, double &eps, double &x0, double x[NMAX+1], int n, int &succes, double RodTabel[NMAX], double &rod){
	char svarNR;
	double test;
	for(int i=0; i<5; i++){
		n = 4;
		x[0] = -10;
		x[1] = 35;
		x[2] = -50;
		x[3] = 24;
		ArrangerKoef(x,n);

		test = 1;

		cout << endl << "Runde nr. " << i+1 << endl;

		cout << "\nUdregner rødder vha. Newton Raphson-metoden: ";
		svarNR = 'j';
		int t = 0;
		if (svarNR=='j' or svarNR=='J'){
			x0 = 0.1;
			cout << endl << "Startpunkt, x0: " << x0;
			eps = 0.001;
			cout << endl << "Epsilon, eps: " << eps;
			N = 1000;
			cout << endl << "Max antal iterationer, N: " << N;
			do {
				NewRapHorner(x0,eps,N,n,x,rod,succes);

				RodTabel[t] = rod;
				t = t+1;
				if (succes == 1){
					HornerDivPol(n,x,rod);
					UdskrivPolynomium(x,n);
					cout << endl << "Sandt resultat: " << setprecision(24) << test;
					cout << endl << "Beregnet resultat: " << setprecision(24) << rod;
					cout << endl << "Præcision af resultat: " << setprecision(24) << test-rod << endl;
					if (n==1){
						rod = -(x[1]/x[0]);
						cout << endl << endl << "Den sidste rod er fundet og vises her: rod = " << rod << endl;
						svarNR = 'n';
						RodTabel[t] = rod;
						cout << endl << "Præcision af resultat: " << setw(10) << test-rod << endl;
					}
					else {
						svarNR = 'j';
						test += 1;
					}
				}
				else {
					cout << "\nNR-iterationen mislykkedes.";
					cout << "\nSe en tabel igen og kontroller, at den viser fortegnsskift";
					cout << "\nEfter konstatering af fortegnsskift: Indtast nye NR-data";
					svarNR = 'n';
				}
			}
			while (svarNR == 'j' or svarNR == 'J');
		}
		cout << endl << "-----------------------------------------------------" << endl;
	}
}




