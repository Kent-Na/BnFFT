/*
 *  BnFFT.h
 *
 *  Created by Kent Nakajima on 09/10/28.
 *
 */

//#include <x86intrin.h>
namespace bn{

//class alignas(128/8) Complex{
class Complex{
public:
	double re;
	double im;
	
	Complex(){};
    //Complex()=default;
	Complex(double reValue);
	Complex(double reValue,double imValue);
	
	inline Complex conjugate(){
		return Complex(re,-im);
	}
	inline void log(){
		printf("%f + i %f\n",re,im);
	}
};

Complex operator + (const Complex &R);
Complex operator - (const Complex &R);

Complex operator + (const Complex &L,const Complex &R);
Complex operator - (const Complex &L,const Complex &R);
Complex operator * (const Complex &L,const Complex &R);

Complex operator * (const Complex &L,const double &R);
Complex operator * (const double &L,const Complex &R);

Complex operator +=	(Complex &L,const Complex &R);
Complex operator -=	(Complex &L,const Complex &R);

inline Complex::Complex(double reValue){
	re=reValue;
	im=0.0;
}
inline Complex::Complex(double reValue,double imValue){
	re=reValue;
	im=imValue;
}

inline Complex operator -	(const Complex &R){
	return Complex(-R.re, -R.im);
}
inline Complex operator +	(const Complex &R){
	return R;
}

inline Complex operator +	(const Complex &L,const Complex &R){
    //__m128d tmp = _mm_add_pd(*(__m128d*)(&L), *(__m128d*)(&R));
    //return *(Complex*)&tmp;
	return Complex(L.re+R.re, L.im+R.im);
}
inline Complex operator -	(const Complex &L,const Complex &R){
    //__m128d tmp = _mm_sub_pd(*(__m128d*)(&L), *(__m128d*)(&R));
    //return *(Complex*)&tmp;
	return Complex(L.re-R.re, L.im-R.im);
}
inline Complex operator *	(const Complex &L,const Complex &R){
    //__m128d tmp_l = *(__m128d*)(&L);
    //__m128d tmp_r = *(__m128d*)(&R);
    //__m128d tmp_l_re = _mm_shuffle_pd(tmp_l, tmp_l, 0);
    //__m128d tmp_l_im = _mm_shuffle_pd(tmp_l, tmp_l, 3);
    //__m128d tmp_r_inv = _mm_shuffle_pd(tmp_l, tmp_l, 2);
    //__m128d tmp_0 = _mm_mul_pd(tmp_l_re, tmp_r);
    //__m128d tmp_1 = _mm_mul_pd(tmp_l_im, tmp_r_inv);
    //__m128d tmp_end = _mm_addsub_pd(tmp_0, tmp_1);
    //return *(Complex*)&tmp_end;
	return Complex(L.re*R.re-L.im*R.im, L.re*R.im+L.im*R.re);
}

inline Complex operator * (const Complex &L,const double &R){
	return Complex(L.re*R, L.im*R);
}

inline Complex operator * (const double &L,const Complex &R){
	return Complex(L*R.re, L*R.im);
}


inline Complex operator +=	(Complex &L,const Complex &R){
	L.re+=R.re;
	L.im+=R.im;
	return L;
}
inline Complex operator -=	(Complex &L,const Complex &R){
	L.re-=R.re;
	L.im-=R.im;
	return L;
}
inline Complex operator *=	(Complex &L,const Complex &R){
	L = L*R;
	return L;
}

template <int depth, bool Inverse = true>
class ComplexFFT2n{
public:
	static const int length=1<<depth;
	
	Complex table[length];
	ComplexFFT2n(){
		int i;
		//make sin,cos table
        for (i=0;i<length;i++){
            if ( not Inverse){
			   table[i]=Complex(cos(i*2*M_PI/length),sin(i*2*M_PI/length));
			}
            else {
			   table[i]=Complex(cos(i*2*M_PI/length),-sin(i*2*M_PI/length));
            }
		}
	}
	void execute(Complex* value){
		Complex tmp;
		int i,j,k;
		int subCount;
		int offset;
		int idush;
		int ioffset;
		int joffset;
		
		//time(sample) culling
		
		for (k=0;k<depth;k++){
			offset=length>>(k+1);
			subCount=1<<k;
			for (j=0;j<subCount;j++){
				joffset=j*offset*2;
				
				//i=0
				{
					idush=joffset;
					ioffset=idush+offset;
					tmp=value[idush]+value[ioffset];
					value[ioffset]=value[idush]-value[ioffset];
					value[idush]=tmp;
				}
				for (i=1;i<offset;i++){
					idush=i+joffset;
					ioffset=idush+offset;
					tmp=value[idush]+value[ioffset];
					value[ioffset]=(value[idush]-value[ioffset])*table[i<<k];
					value[idush]=tmp;
				}
			}
		}
		
		int n=length;
		i = 0;
		for (j = 1; j < n - 1; j++) {
			for (k = n >> 1; k > (i ^= k); k >>= 1);
			if (j < i) {
				tmp = value[j];
				value[j] = value[i];
				value[i] = tmp;
			}
		}
		
	}
};
template <int depth>
class ComplexFFT2n3{
	
public:
	static const int N1=3;
	static const int N2=(1<<depth);
	static const int length=N1*N2;
	
	Complex table[N2];
	
	ComplexFFT2n3(bool inverse=false){
		int i;
		//make sin,cos table
        for (i=0;i<N2;i++){
            if (inverse){
				table[i]=Complex(cos(i*2*M_PI/N2),-sin(i*2*M_PI/N2));
			}
            else{
				table[i]=Complex(cos(i*2*M_PI/N2),sin(i*2*M_PI/N2));
			}
		}
	}
	
	void execute(Complex* value){
		//3 2^n FFT
		//PrimeFactor
		
		//N1 = 3
		//N2 = 2^n
		//N  = N1 x N2
		
		//t1 = 171 ,t2 = 1
		//t1 N1 + t2 N2 = N +1
		
		int i;
		
		const int t1=(((2-(depth&1))*N2)+1)/3;
		const int t2=1+(depth&1);
		
		//3 sample fft
		{
			int i1,i2,i3;
			Complex tmp1,tmp2,tmp3;
			const Complex W31(0, sqrt(3)*0.5);
			if (t2==1){
				for (i=0;i<N2;i++){
					i1=(N1*t1*i)%length;
					i2=(i1+N2*t2)%length;
					i3=(i2+N2*t2)%length;
					
					tmp1=value[i2]+value[i3];
					tmp2=value[i1]-tmp1*0.5;
					tmp3=(value[i2]-value[i3])*W31;
					
					value[i2]=tmp2+tmp3;
					value[i3]=tmp2-tmp3;
					value[i1]+=tmp1;
				}
			}
			else{
				for (i=0;i<N2;i++){
					i1=(N1*t1*i)%length;
					i2=(i1+N2*t2)%length;
					i3=(i2+N2*t2)%length;
					
					tmp1=value[i2]+value[i3];
					tmp2=value[i1]-tmp1*0.5;
					tmp3=(value[i3]-value[i2])*W31;
					
					value[i2]=tmp2+tmp3;
					value[i3]=tmp2-tmp3;
					value[i1]+=tmp1;
				}
			}
		}	
		//2^n sample fft
		{
			Complex tmp;
			int i,j,k,l;
			int subCount;
			int offset;
			int idush;
			int jdush;
			int ioffset;
			int joffset;
			
			//time(sample) culling
			
			for (l=0;l<N1;l++){	
				for (k=0;k<depth;k++){
					offset=N2>>(k+1);
					subCount=1<<k;
					for (j=0;j<subCount;j++){
						joffset=j*offset*2;
						
						//i=0
						{
							idush=joffset;
							ioffset=idush+offset;
							idush=(idush*N1*t1+N2*t2*l)%length;
							ioffset=(ioffset*N1*t1+N2*t2*l)%length;
							tmp=value[idush]+value[ioffset];
							value[ioffset]=value[idush]-value[ioffset];
							value[idush]=tmp;
						}
						for (i=1;i<offset;i++){
							idush=i+joffset;
							ioffset=idush+offset;
							idush=(idush*N1*t1+N2*t2*l)%length;
							ioffset=(ioffset*N1*t1+N2*t2*l)%length;
							tmp=value[idush]+value[ioffset];
							value[ioffset]=(value[idush]-value[ioffset])*table[((i<<k)*t1)%N2];
							value[idush]=tmp;
						}
					}
				}
			}
			
			int n=N2;
			i = 0;
			for (j = 1; j < n - 1; j++) {
				for (k = n >> 1; k > (i ^= k); k >>= 1);
				if (j < i) {
					for (l=0;l<N1;l++){			
						idush=(i*N1*t1+N2*t2*l)%length;
						jdush=(j*N1*t1+N2*t2*l)%length;
						
						tmp = value[jdush];
						value[jdush] = value[idush];
						value[idush] = tmp;
					}
				}
			}
		}
	}
};

template <int depth, bool Inverse = false, typename baseFFT = ComplexFFT2n<depth-1, Inverse> >
class RealFFT{
	baseFFT cFFT;

public:
	static const int length=1<<depth;	
	Complex table[length];

	RealFFT(){
		int i;
		//make sin,cos table
        for (i=0;i<length;i++){
            if ( not Inverse){
				table[i]=
                    Complex(cos(i*2*M_PI/length),sin(i*2*M_PI/length));
			}
            else {
				table[i]=
                    Complex(cos(i*2*M_PI/length),-sin(i*2*M_PI/length));
			}
		}
	}
	void execute(double* value){
		Complex* cValue = (Complex*)value;
		cFFT.execute(cValue);
		const int halfLength = length/2;
		
		for (int i = 1; i<halfLength/2; i++){
			Complex c0 = cValue[i];
			Complex c1 = cValue[halfLength-i];
			
			cValue[i]=c0-(0.5+Complex(0, 0.5)*table[i])*(c0-c1.conjugate());
			cValue[halfLength-i]=c1-(0.5+Complex(0, 0.5)*table[halfLength-i])*(c1-c0.conjugate());
		}
		{
			Complex c0 = cValue[0];
			value[0]=c0.re+c0.im;
			value[1]=c0.re-c0.im;
		}
	}
	
};

template <typename baseFFT>
class Complex2DFFT{
	baseFFT FFT;
public:
	void execute(Complex* value,unsigned int width){
		int i,j;
		for (i=0;i<FFT.length;i++){
			FFT.execute(value+i*width);
		}
		
		//transpose
		Complex tmp;
		for (i=0;i<FFT.length;i++){
			for (j=0;j<i;j++){
				tmp=value[i+j*width];
				value[i+j*width]=value[j+i*width];
				value[j+i*width]=tmp;
			}
		}
		
		for (i=0;i<FFT.length;i++){
			FFT.execute(value+i*width);
		}
		
		//inv transpose(option)
		if (0){
			Complex tmp;
			for (i=0;i<FFT.length;i++){
				for (j=0;j<i;j++){
					tmp=value[i+j*width];
					value[i+j*width]=value[j+i*width];
					value[j+i*width]=tmp;
				}
			}
		}
	}
};

}// End of namespace
