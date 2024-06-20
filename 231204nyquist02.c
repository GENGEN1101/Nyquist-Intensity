//スプリットステップフーリエ法
//構造体、複素配列complex型、数値計算ライブラリFFTW3利用
//伝送路計算プログラム
//150910 パルス形状をナイキストパルスに変更 ガウスやsech型は"140429ssfｍv3.c"を参照
//211115 変調器と2本のファイバを通る
//211217 gain追加, Capectrumf追加(横軸周波数で出力される), Modulatorの変数B 修正
//220608 specout01,02のみの出力, 損失・Ds_1・beta3_1を考慮
//220704 コアダンプ発生はスペクトルの計算s領域を超えているから
//220704 grid.lambdaとlambda_Nとlambda0_1の値を同じにするとコアダンプがなおる可能性あり
//220908 損失・Ds_1・beta3_1を考慮しない	//beta2を変更して比較
//221029 DCF考慮(beta2の値は手入力)
//221204 コアダンプ修正済み

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<complex.h>
#include<fftw3.h>
#include"funcv3.h"

#define PI 3.1415926535897932384626
#define CV 2.99792458
#define N 1024 		/* point number of window ＦＦＴのポイント数、計算のため(2^n)の数値とすること*/
#define NUM_FIBER 100		/* スプリットステップファイバ計算分割数 計算精度に依存するので小さすぎる値にしないこと*/
#define WIDTH_WINDOW 100.0	/* ウィンドウ（計算領域）幅[ps] パルスの裾が両端にかからないようにすること */

int Modulator(struct OPTSIG* optsig, struct GRID* grid, double A, double B, double To, double Ep1);
int Gain(struct OPTSIG* optsig, struct GRID* grid, double gain);
int Cspectrumf(FILE *fp, struct OPTSIG* optsig, struct GRID* grid);

/* メイン関数 */
int main()
{
//	int i;
	//各構造体の準備（定義はfuncv3.h参照)
	struct GRID grid;		//時間、周波数グリッド用構造体
	struct FIBER fiber_1, fiber_2;		//ファイバパラメータ用構造体
	struct OPTSIG optsig;	//光信号用構造体
	struct OPTSIG o;
	
	double Z;

	char tfile[120],wfile[120];
	FILE *CSPEC1,*CWAVE1;
	char cspecfile1[120],cwavefile1[120];
	
	// ========== Simulation parameter 構造体GRIDの初期設定 ==========
	grid.lambda = 1.550;			// Central wavelength in simulation[micro meter]、数値シミュレーションの中心波長
	grid.window = WIDTH_WINDOW;		// time window
	grid.n = N;						// FFT size
	grid.dt =  WIDTH_WINDOW/N;		// time step size
	// =========================================
	
	// ========== Nyquistpulse parameter  ==========
	double lambda_N = grid.lambda;				/* Wavelength [micro meter] */
	double Pp_N = 1.514;					/* ナイキストパルスのピーク電力 [W]平均電力で指定の場合無視される */
	double timeslot_N = WIDTH_WINDOW/4.0;			/* timeslot [ps] */
	double alpha_N = 0.75;				/* roll-off factor */
	double delaytime_N =0.0;				/* 郡遅延 */
	// =============================================
	
	// ========== Signal parameter  ==========
	double lambda1 = grid.lambda;			/* Wavelength [micro meter] grid.lambdaと同じ値にすることでコアダンプ発生を防ぐ*/
	double Pp1 = 0.670;					/* Power [W] */	
	//double Tfwhm1 = 1.0;				/* FWHM time width [ps] <=0 なら連続光 */
	// =============================================
	double Ep1 = sqrt(Pp1);

	//==========Modulater=============
	double m = 0.80;	//変調度
	double A = m*PI/4.0;     //振幅
	double B = A-PI/2.0;     //直流成分
	double To= WIDTH_WINDOW;      //変調器の周期(WIDTH_WINDOWと同じにすればうまく表示される)
	//==============================

	//double LD=1.0;
	//double K_1 = 66.0;            //LD/LNL=K^2
	//double LNL_1 = 1.0 / (K_1*K_1);
	
	/* ==========Fiber parameter1========== */
	double Fiberloss_1 = 0.0; 				/* Fiber Loss [dB/km] */
	//double lambda0_1 = 1.567749;				/* zero dispersion wavelength [micro meter] Dを直接指定する場合は無視 */
	//double Ds_1 = 0.02077;						/* dispersion slope [ps/km nm^2] */
	//double D_1 = Ds_1*(grid.lambda-lambda0_1)*1000;				/*分散スロープを考慮*/
	//double D_1 = -0.391856;						/*Dispersion parameter [ps/ km nm]*/
	/*		2*PI*CV/(lambda1*lambda1) = 0.784580892 [1/(ps nm)] */
	double gamma_1 = 12;					/* nonlinear parameter [1/(W km)] */
	double length_1 = 0.403;				/* Fiber length [km] */
	/* ====================================== */
		
	//double beta2_1 = 0.391856 /0.784580892   ;					/* (Disp para [ps/ km nm])/0.784580892		 [ps^2/ km] */
	double beta2_1 = 0.496   ;
	double beta3_1 = 0.0; /* 100.0*pow(grid.lambda*grid.lambda / (2.0*PI*CV), 2)*(Ds_2 + D_2 / (500.0*grid.lambda));	[ps^3/km]  */
	//double beta2_1 = -10.0*D_1*grid.lambda*grid.lambda/(2.0*PI*CV);							/* [ps^2/ km] */
	//double beta3_1 =100.0*pow(grid.lambda*grid.lambda/(2.0*PI*CV),2)*(Ds_1+D_1/(500.0*grid.lambda));	/* [ps^3/km]  */
	double alpha_1  = 0.115129254*Fiberloss_1;
	double dz_1 =length_1/NUM_FIBER;
	
	/* ==========Fiber parameter2==========（分散補償用） */
	double Fiberloss_2 = 0.0; 				/* Fiber Loss [dB/km] */
	//double lambda0_2 = 1.550;				/* zero dispersion wavelength [micro meter] Dを直接指定する場合は無視  */
	//double Ds_2 = 0.0;						/* dispersion slope [ps/km nm^2] */
	//double D_2 = 0.006951808;				/*Ds*( grid.lambda - lambda0)*1000 Dispersion parameter [ps/ km nm] Ds*(lambda-lambda0)*1000*/
	double gamma_2 = 0.0;						/* nonlinear parameter [1/(W km)] */
	double length_2 = 9.77 ;						/* Fiber length [km] */
	/* ====================================== */
	
	double beta2_2 = -16.0   ;					/* (Disp para [ps/ km nm])/0.784580892		 [ps^2/ km] */
	//double beta2_2 = -10.0*D_2*grid.lambda*grid.lambda/(2.0*PI*CV);							/* [ps^2/ km] */
	double beta3_2 = 0.0; /* 100.0*pow(grid.lambda*grid.lambda / (2.0*PI*CV), 2)*(Ds_2 + D_2 / (500.0*grid.lambda));	[ps^3/km]  */
	//double beta3_2 = 100.0*pow(grid.lambda*grid.lambda/(2.0*PI*CV),2)*(Ds_2+D_2/(500.0*grid.lambda));	/* [ps^3/km]  */
	double alpha_2 = 0.115129254*Fiberloss_2;
	double dz_2 =length_2/NUM_FIBER;	
	
	
	// gridに設定したパラメータを用い光波形用構造体optsigの設定、波形メモリ確保、FFTW3のplanを設定
	Initfftw(&optsig, &grid);
	Initfftw(&o, &grid);
	// Ssfm用ファイバパラメータ構造体設定,beta計算用メモリ確保
	Initfiber(&fiber_1, &grid);
	Fiberparam(&fiber_1, &grid, beta2_1, beta3_1, alpha_1, gamma_1, dz_1, length_1, NUM_FIBER);
	Initfiber(&fiber_2, &grid);
	Fiberparam(&fiber_2, &grid, beta2_2, beta3_2, alpha_2, gamma_2, dz_2, length_2, NUM_FIBER);	
	

	/* 伝送路中の波形、スペクトルのファイル名設定 */
	//sprintf(tfile,"tzq_c1.csv");
	//sprintf(wfile,"wzq_c1.csv");
	
	/* ========== initial waveform ========== */
	/* CW光源(逆フーリエ変換) */
	Cwsource(&optsig,&grid,lambda1,Ep1, 1);	//アドレス渡しのため&をつけている
	//Cwsource(&o, &grid, lambda1, Ep1, 1);	//アドレス渡しのため&をつけている


	//coswave(&optsig, &grid);
	Modulator(&optsig, &grid, A, B, To, Ep1);
	
	//増幅器
	double G = 1 / (cos(A + B)*cos(A + B));
	
	if (m > 1){
		G = 1 / ((cos(A + B)*cos(A + B)+1)/2);
	}
	
	Gain(&optsig, &grid, G);

	
	/* ==========目標波形(ナイキストパルス)========== */
	RootNyquistpulse(&o, &grid, lambda_N, timeslot_N, alpha_N, Pp_N, delaytime_N);//ルートナイキストパルス（信号光）
	//Nyquistpulse(&o, &grid, lambda_N, timeslot_N, alpha_N, Pp_N, delaytime_N);	//ナイキストパルス（信号光）

	
	
	
	/* 伝送前波形、スペクトルを出力 */
	/*
	sprintf(cwavefile1,"cwaveinc1.csv");
	CWAVE1 = fopen(cwavefile1,"w");
	fprintf(CWAVE1, "# Time [ps], P[W] Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Cwaveform(CWAVE1,&optsig, &grid);
	fclose(CWAVE1);
	
	sprintf(cspecfile1, "cspecinc1.csv");
	CSPEC1 = fopen(cspecfile1, "w");
	fprintf(CSPEC1, "# Lambda [nm],f [THz], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Delay(&optsig, &grid, lambda1, -1.0*WIDTH_WINDOW / 2.0);
	Cspectrumf(CSPEC1, &optsig, &grid);
	Delay(&optsig, &grid, lambda1, WIDTH_WINDOW / 2.0);
	fclose(CSPEC1);
	*/
	/*
	sprintf(cspecfile1,"cspecinc1.csv");
	CSPEC1 = fopen(cspecfile1,"w");
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1点あたりの電力、PSDが欲しいときは単位は変換が必要
	Delay(&optsig,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//スペクトルを見やすくするため群遅延時間を-1.0*WIDTH_WINDOW/2.0与える。（位相をフラットにする）
	Cspectrum(CSPEC1,&optsig, &grid);
	Delay(&optsig,&grid,lambda1,WIDTH_WINDOW/2.0);			//群遅延時間を1.0*WIDTH_WINDOW/2.0与え、元に戻す。
	fclose(CSPEC1);
	*/
	
	/* スプリットステップフーリエ法 */
	Z = 0.0;	//初期ファイバ長リセット 	
	//Ssfm_fo(&optsig, &grid, &fiber_1, &Z, tfile, wfile);		//伝搬中波形ファイル出力有 １本目
	Ssfm_nfo(&optsig, &grid, &fiber_1, &Z);						//伝搬中波形ファイル出力無し 1本目
	
	
	//Ssfm_fo(&optsig, &grid, &fiber_2, &Z, tfile, wfile);		//伝搬中波形ファイル出力有 ２本目
	//Ssfm_nfo(&optsig, &grid, &fiber, &Z);					//伝搬中波形ファイル出力無の場合の例

	/* 1段目伝送後波形、スペクトルを出力 */
	sprintf(cwavefile1,"cwaveoutc01_beta2=%0.3f_m=%0.1f.csv", beta2_1, m);
	CWAVE1 = fopen(cwavefile1,"w");
	fprintf(CWAVE1, "# Time [ps], P[W] Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Cwaveform(CWAVE1,&optsig, &grid);
	fclose(CWAVE1);
	
	sprintf(cspecfile1, "cspecoutc01_beta2=%0.3f_m=%0.1f.csv", beta2_1, m);
	CSPEC1 = fopen(cspecfile1, "w");
	fprintf(CSPEC1, "# Lambda [nm],f [THz], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Delay(&optsig, &grid, lambda1, -1.0*WIDTH_WINDOW / 2.0);
	Cspectrumf(CSPEC1, &optsig, &grid);
	Delay(&optsig, &grid, lambda1, WIDTH_WINDOW / 2.0);
	fclose(CSPEC1);
	
	/*
	sprintf(cspecfile1,"cspecoutc01.csv");
	CSPEC1 = fopen(cspecfile1,"w");
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1点あたりの電力、PSDが欲しいときは単位は変換が必要
	Delay(&optsig,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//スペクトルを見やすくするため群遅延時間を-1.0*WIDTH_WINDOW/2.0与える。（位相をフラットにする）
	Cspectrum(CSPEC1,&optsig, &grid);
	Delay(&optsig,&grid,lambda1,WIDTH_WINDOW/2.0);			//群遅延時間をWIDTH_WINDOW/2.0与え、元に戻す。
	fclose(CSPEC1);
	*/
	
	/* スプリットステップフーリエ法 */
	//Ssfm_fo(&optsig, &grid, &fiber_2, &Z, tfile, wfile);
	Ssfm_nfo(&optsig, &grid, &fiber_2, &Z);
	
	/* 2段目伝送後波形、スペクトルを出力 */
	sprintf(cwavefile1,"DCFwaveoutc_beta2=%0.3f_m=%0.1f.csv", beta2_1, m);
	CWAVE1 = fopen(cwavefile1,"w");
	fprintf(CWAVE1, "# Time [ps], P[W] Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Cwaveform(CWAVE1,&optsig, &grid);
	fclose(CWAVE1);
	
	sprintf(cspecfile1, "DCFspecoutc_beta2=%0.3f_m=%0.1f.csv", beta2_1, m);
	CSPEC1 = fopen(cspecfile1, "w");
	fprintf(CSPEC1, "# Lambda [nm],f [THz], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Delay(&optsig, &grid, lambda1, -1.0*WIDTH_WINDOW / 2.0);
	Cspectrumf(CSPEC1, &optsig, &grid);
	Delay(&optsig, &grid, lambda1, WIDTH_WINDOW / 2.0);
	fclose(CSPEC1);
	
	/*
	sprintf(cspecfile1,"cspecoutc02.csv");
	CSPEC1 = fopen(cspecfile1,"w");
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1点あたりの電力、PSDが欲しいときは単位は変換が必要
	Delay(&optsig,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//スペクトルを見やすくするため群遅延時間を-1.0*WIDTH_WINDOW/2.0与える。（位相をフラットにする）
	Cspectrum(CSPEC1,&optsig, &grid);
	Delay(&optsig,&grid,lambda1,WIDTH_WINDOW/2.0);			//群遅延時間をWIDTH_WINDOW/2.0与え、元に戻す。
	fclose(CSPEC1);
	*/
	
	/* 目標出力(ナイキストパルス) */ 
	sprintf(cwavefile1,"nyq_wave.csv");
	CWAVE1 = fopen(cwavefile1,"w");
	fprintf(CWAVE1, "# Time [ps], P[W] Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Cwaveform(CWAVE1,&o, &grid);
	fclose(CWAVE1);
	
	sprintf(cspecfile1, "nyq_spec.csv");
	CSPEC1 = fopen(cspecfile1, "w");
	fprintf(CSPEC1, "# Lambda [nm],f [THz], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");
	Delay(&o, &grid, lambda1, -1.0*WIDTH_WINDOW / 2.0);
	Cspectrumf(CSPEC1, &o, &grid);
	Delay(&o, &grid, lambda1, WIDTH_WINDOW / 2.0);
	fclose(CSPEC1);
	
	/*
	sprintf(cspecfile1,"nyq_spec.csv");
	CSPEC1 = fopen(cspecfile1,"w");
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1点あたりの電力、PSDが欲しいときは単位は変換が必要
	Delay(&o,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//スペクトルを見やすくするため群遅延時間を-1.0*WIDTH_WINDOW/2.0与える。（位相をフラットにする）
	Cspectrum(CSPEC1,&o, &grid);
	Delay(&o,&grid,lambda1,WIDTH_WINDOW/2.0);			//群遅延時間をWIDTH_WINDOW/2.0与え、元に戻す。
	fclose(CSPEC1);
	*/
	
	Closefftw(&optsig);
	Closefftw(&o);
	Closefiber(&fiber_1);
	Closefiber(&fiber_2);
	return(0);
}

int Modulator(struct OPTSIG* optsig, struct GRID* grid, double A, double B, double To,double Ep1)
{
	int j, n;
	double f = 1.0/To;//周波数
	double W, dT, T;
	complex *q;


	q = optsig->q;
	W = grid->window;
	n = grid->n;
	dT = grid->dt;

	T = W / n;


	for (j = 0; j < n; j++) {
		q[j] = Ep1 * cos(A*cos(2 * PI*f*(j - n / 2)*T) + B);
	}
	return(0);
}

//*増幅器　*/
int Gain(struct OPTSIG* optsig, struct GRID* grid, double gain)
{
	int i;
	double egain = sqrt(gain);

	for (i = 0; i < grid->n; i++) {
		optsig->q[i] = optsig->q[i] * egain;
	}

	return(0);
}

int Cspectrumf(FILE *fp, struct OPTSIG* optsig, struct GRID* grid)
{
	int k;
	double f, f0, lambda, window;//w,
	int n;

	n = grid->n;
	window = grid->window;

	struct OPTSIG optsig2;
	Initfftw(&optsig2, grid);


	memcpy(optsig2.q, optsig->q, sizeof(fftw_complex)*n);
	fftw_execute(optsig2.pfft);

	f0 = 100.0*CV / grid->lambda; /* THz */
	

	for (k = n / 2; k <= n - 1; k = k + 1)
	{
		//		w = 2.0*PI*(k- n)/window;
		f = (k - n) / window + f0;	/* THz */
		lambda = 100.0*CV / f; 	/* micro m */
		fprintf(fp, "%e, %e, %e, %e, %e, %e\n", 1000 * lambda, f-f0, creal(optsig2.q[k] * conj(optsig2.q[k]) / n / n), creal(optsig2.q[k]), cimag(optsig2.q[k]), carg(optsig2.q[k]));
		/* nm,  W, Real, Imaginaly, phase */
	}
	for (k = 0; k <= n / 2 - 1; k = k + 1)
	{
		//		w = 2.0*PI*k/window;
		f = k / window + f0; /* THz */
		lambda = 100.0*CV / f; /* micro m */
		fprintf(fp, "%e, %e, %e, %e, %e,%e\n", 1000 * lambda, f - f0, creal(optsig2.q[k] * conj(optsig2.q[k]) / n / n), creal(optsig2.q[k]), cimag(optsig2.q[k]), carg(optsig2.q[k]));
		/* nm, W, Real, Imaginaly, phase */
	}
	fprintf(fp, "\n\n");
	Closefftw(&optsig2);
	return(0);
}

