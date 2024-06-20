//�X�v���b�g�X�e�b�v�t�[���G�@
//�\���́A���f�z��complex�^�A���l�v�Z���C�u����FFTW3���p
//�`���H�v�Z�v���O����
//150910 �p���X�`����i�C�L�X�g�p���X�ɕύX �K�E�X��sech�^��"140429ssf��v3.c"���Q��
//211115 �ϒ����2�{�̃t�@�C�o��ʂ�
//211217 gain�ǉ�, Capectrumf�ǉ�(�������g���ŏo�͂����), Modulator�̕ϐ�B �C��
//220608 specout01,02�݂̂̏o��, �����EDs_1�Ebeta3_1���l��
//220704 �R�A�_���v�����̓X�y�N�g���̌v�Zs�̈�𒴂��Ă��邩��
//220704 grid.lambda��lambda_N��lambda0_1�̒l�𓯂��ɂ���ƃR�A�_���v���Ȃ���\������
//220908 �����EDs_1�Ebeta3_1���l�����Ȃ�	//beta2��ύX���Ĕ�r
//221029 DCF�l��(beta2�̒l�͎����)
//221204 �R�A�_���v�C���ς�

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<complex.h>
#include<fftw3.h>
#include"funcv3.h"

#define PI 3.1415926535897932384626
#define CV 2.99792458
#define N 1024 		/* point number of window �e�e�s�̃|�C���g���A�v�Z�̂���(2^n)�̐��l�Ƃ��邱��*/
#define NUM_FIBER 100		/* �X�v���b�g�X�e�b�v�t�@�C�o�v�Z������ �v�Z���x�Ɉˑ�����̂ŏ���������l�ɂ��Ȃ�����*/
#define WIDTH_WINDOW 100.0	/* �E�B���h�E�i�v�Z�̈�j��[ps] �p���X�̐������[�ɂ�����Ȃ��悤�ɂ��邱�� */

int Modulator(struct OPTSIG* optsig, struct GRID* grid, double A, double B, double To, double Ep1);
int Gain(struct OPTSIG* optsig, struct GRID* grid, double gain);
int Cspectrumf(FILE *fp, struct OPTSIG* optsig, struct GRID* grid);

/* ���C���֐� */
int main()
{
//	int i;
	//�e�\���̂̏����i��`��funcv3.h�Q��)
	struct GRID grid;		//���ԁA���g���O���b�h�p�\����
	struct FIBER fiber_1, fiber_2;		//�t�@�C�o�p�����[�^�p�\����
	struct OPTSIG optsig;	//���M���p�\����
	struct OPTSIG o;
	
	double Z;

	char tfile[120],wfile[120];
	FILE *CSPEC1,*CWAVE1;
	char cspecfile1[120],cwavefile1[120];
	
	// ========== Simulation parameter �\����GRID�̏����ݒ� ==========
	grid.lambda = 1.550;			// Central wavelength in simulation[micro meter]�A���l�V�~�����[�V�����̒��S�g��
	grid.window = WIDTH_WINDOW;		// time window
	grid.n = N;						// FFT size
	grid.dt =  WIDTH_WINDOW/N;		// time step size
	// =========================================
	
	// ========== Nyquistpulse parameter  ==========
	double lambda_N = grid.lambda;				/* Wavelength [micro meter] */
	double Pp_N = 1.514;					/* �i�C�L�X�g�p���X�̃s�[�N�d�� [W]���ϓd�͂Ŏw��̏ꍇ��������� */
	double timeslot_N = WIDTH_WINDOW/4.0;			/* timeslot [ps] */
	double alpha_N = 0.75;				/* roll-off factor */
	double delaytime_N =0.0;				/* �S�x�� */
	// =============================================
	
	// ========== Signal parameter  ==========
	double lambda1 = grid.lambda;			/* Wavelength [micro meter] grid.lambda�Ɠ����l�ɂ��邱�ƂŃR�A�_���v������h��*/
	double Pp1 = 0.670;					/* Power [W] */	
	//double Tfwhm1 = 1.0;				/* FWHM time width [ps] <=0 �Ȃ�A���� */
	// =============================================
	double Ep1 = sqrt(Pp1);

	//==========Modulater=============
	double m = 0.80;	//�ϒ��x
	double A = m*PI/4.0;     //�U��
	double B = A-PI/2.0;     //��������
	double To= WIDTH_WINDOW;      //�ϒ���̎���(WIDTH_WINDOW�Ɠ����ɂ���΂��܂��\�������)
	//==============================

	//double LD=1.0;
	//double K_1 = 66.0;            //LD/LNL=K^2
	//double LNL_1 = 1.0 / (K_1*K_1);
	
	/* ==========Fiber parameter1========== */
	double Fiberloss_1 = 0.0; 				/* Fiber Loss [dB/km] */
	//double lambda0_1 = 1.567749;				/* zero dispersion wavelength [micro meter] D�𒼐ڎw�肷��ꍇ�͖��� */
	//double Ds_1 = 0.02077;						/* dispersion slope [ps/km nm^2] */
	//double D_1 = Ds_1*(grid.lambda-lambda0_1)*1000;				/*���U�X���[�v���l��*/
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
	
	/* ==========Fiber parameter2==========�i���U�⏞�p�j */
	double Fiberloss_2 = 0.0; 				/* Fiber Loss [dB/km] */
	//double lambda0_2 = 1.550;				/* zero dispersion wavelength [micro meter] D�𒼐ڎw�肷��ꍇ�͖���  */
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
	
	
	// grid�ɐݒ肵���p�����[�^��p�����g�`�p�\����optsig�̐ݒ�A�g�`�������m�ہAFFTW3��plan��ݒ�
	Initfftw(&optsig, &grid);
	Initfftw(&o, &grid);
	// Ssfm�p�t�@�C�o�p�����[�^�\���̐ݒ�,beta�v�Z�p�������m��
	Initfiber(&fiber_1, &grid);
	Fiberparam(&fiber_1, &grid, beta2_1, beta3_1, alpha_1, gamma_1, dz_1, length_1, NUM_FIBER);
	Initfiber(&fiber_2, &grid);
	Fiberparam(&fiber_2, &grid, beta2_2, beta3_2, alpha_2, gamma_2, dz_2, length_2, NUM_FIBER);	
	

	/* �`���H���̔g�`�A�X�y�N�g���̃t�@�C�����ݒ� */
	//sprintf(tfile,"tzq_c1.csv");
	//sprintf(wfile,"wzq_c1.csv");
	
	/* ========== initial waveform ========== */
	/* CW����(�t�t�[���G�ϊ�) */
	Cwsource(&optsig,&grid,lambda1,Ep1, 1);	//�A�h���X�n���̂���&�����Ă���
	//Cwsource(&o, &grid, lambda1, Ep1, 1);	//�A�h���X�n���̂���&�����Ă���


	//coswave(&optsig, &grid);
	Modulator(&optsig, &grid, A, B, To, Ep1);
	
	//������
	double G = 1 / (cos(A + B)*cos(A + B));
	
	if (m > 1){
		G = 1 / ((cos(A + B)*cos(A + B)+1)/2);
	}
	
	Gain(&optsig, &grid, G);

	
	/* ==========�ڕW�g�`(�i�C�L�X�g�p���X)========== */
	RootNyquistpulse(&o, &grid, lambda_N, timeslot_N, alpha_N, Pp_N, delaytime_N);//���[�g�i�C�L�X�g�p���X�i�M�����j
	//Nyquistpulse(&o, &grid, lambda_N, timeslot_N, alpha_N, Pp_N, delaytime_N);	//�i�C�L�X�g�p���X�i�M�����j

	
	
	
	/* �`���O�g�`�A�X�y�N�g�����o�� */
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
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1�_������̓d�́APSD���~�����Ƃ��͒P�ʂ͕ϊ����K�v
	Delay(&optsig,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//�X�y�N�g�������₷�����邽�ߌQ�x�����Ԃ�-1.0*WIDTH_WINDOW/2.0�^����B�i�ʑ����t���b�g�ɂ���j
	Cspectrum(CSPEC1,&optsig, &grid);
	Delay(&optsig,&grid,lambda1,WIDTH_WINDOW/2.0);			//�Q�x�����Ԃ�1.0*WIDTH_WINDOW/2.0�^���A���ɖ߂��B
	fclose(CSPEC1);
	*/
	
	/* �X�v���b�g�X�e�b�v�t�[���G�@ */
	Z = 0.0;	//�����t�@�C�o�����Z�b�g 	
	//Ssfm_fo(&optsig, &grid, &fiber_1, &Z, tfile, wfile);		//�`�����g�`�t�@�C���o�͗L �P�{��
	Ssfm_nfo(&optsig, &grid, &fiber_1, &Z);						//�`�����g�`�t�@�C���o�͖��� 1�{��
	
	
	//Ssfm_fo(&optsig, &grid, &fiber_2, &Z, tfile, wfile);		//�`�����g�`�t�@�C���o�͗L �Q�{��
	//Ssfm_nfo(&optsig, &grid, &fiber, &Z);					//�`�����g�`�t�@�C���o�͖��̏ꍇ�̗�

	/* 1�i�ړ`����g�`�A�X�y�N�g�����o�� */
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
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1�_������̓d�́APSD���~�����Ƃ��͒P�ʂ͕ϊ����K�v
	Delay(&optsig,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//�X�y�N�g�������₷�����邽�ߌQ�x�����Ԃ�-1.0*WIDTH_WINDOW/2.0�^����B�i�ʑ����t���b�g�ɂ���j
	Cspectrum(CSPEC1,&optsig, &grid);
	Delay(&optsig,&grid,lambda1,WIDTH_WINDOW/2.0);			//�Q�x�����Ԃ�WIDTH_WINDOW/2.0�^���A���ɖ߂��B
	fclose(CSPEC1);
	*/
	
	/* �X�v���b�g�X�e�b�v�t�[���G�@ */
	//Ssfm_fo(&optsig, &grid, &fiber_2, &Z, tfile, wfile);
	Ssfm_nfo(&optsig, &grid, &fiber_2, &Z);
	
	/* 2�i�ړ`����g�`�A�X�y�N�g�����o�� */
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
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1�_������̓d�́APSD���~�����Ƃ��͒P�ʂ͕ϊ����K�v
	Delay(&optsig,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//�X�y�N�g�������₷�����邽�ߌQ�x�����Ԃ�-1.0*WIDTH_WINDOW/2.0�^����B�i�ʑ����t���b�g�ɂ���j
	Cspectrum(CSPEC1,&optsig, &grid);
	Delay(&optsig,&grid,lambda1,WIDTH_WINDOW/2.0);			//�Q�x�����Ԃ�WIDTH_WINDOW/2.0�^���A���ɖ߂��B
	fclose(CSPEC1);
	*/
	
	/* �ڕW�o��(�i�C�L�X�g�p���X) */ 
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
	fprintf(CSPEC1, "# Lambda [nm], P[W], Real[W^0.5], Imaginaly[W^0.5], Phase[rad], #\n");//1�_������̓d�́APSD���~�����Ƃ��͒P�ʂ͕ϊ����K�v
	Delay(&o,&grid,lambda1,-1.0*WIDTH_WINDOW/2.0);		//�X�y�N�g�������₷�����邽�ߌQ�x�����Ԃ�-1.0*WIDTH_WINDOW/2.0�^����B�i�ʑ����t���b�g�ɂ���j
	Cspectrum(CSPEC1,&o, &grid);
	Delay(&o,&grid,lambda1,WIDTH_WINDOW/2.0);			//�Q�x�����Ԃ�WIDTH_WINDOW/2.0�^���A���ɖ߂��B
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
	double f = 1.0/To;//���g��
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

//*������@*/
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

