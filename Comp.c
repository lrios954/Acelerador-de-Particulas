#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define q 1.602
#define m 1.672
#define B 0.03
#define d 6
#define R 800
#define a 20

//a viene definido por el campo electrico, la carga y la masa
int main() {

//Archivo donde se guardaran los datos de la trayectoria
FILE *export;
  export = fopen("trayectoria.dat", "w");
  double h;
  int n_points;
  h = 1;
  n_points = 200;
  
  
  double *t;
  double *Vel;
  double *Pos;


 double *theta;
 double *Costheta;
 double *Sentheta;

double *Posx;
double *Posy;  
  
  Vel = malloc(n_points*sizeof(double));
  Pos = malloc(n_points*sizeof(double));
  t = malloc(n_points*sizeof(double));
  theta = malloc(n_points*sizeof(double));
  Costheta = malloc(n_points*sizeof(double));
  Sentheta = malloc(n_points*sizeof(double));
  Posx= malloc(n_points*sizeof(double));
  Posy= malloc(n_points*sizeof(double));


  
  if (!Vel || !Pos || !t){
    printf("Error en memoria");
    exit(1);
  }
  
  
//condiciones iniciales del problema

  Pos[0] = 1.0;
  t[0] = 0.0;
  theta[0]=0.0;
if(Pos[0]>0){
Costheta[0]=0.1;
}
else if(Pos[0]<0){
Costheta[0]=-0.1;
}
  
Sentheta[0]=sin(theta[0]);
if(Costheta[0]>0){
	
	Vel[0]=sqrt(2*d*a);
			
}
else if(Costheta[0]<0){
			
	Vel[0]=-sqrt(2*d*(a));
			
}

  
  fprintf(export,"%f %f %f %f %f \n", t[0],Vel[0],Pos[0],Costheta[0], Sentheta[0]);

  int i;
  
  for (i = 1; i < n_points; i ++)
    {
      rungekutta(i,Vel,Pos,t,h);
      theta[i]=10*i/(57.2956);
      Costheta[i]=cos(theta[i]);
      Sentheta[i]=sin(theta[i]);

      fprintf(export,"%f %f %f %f %f \n",t[i],Vel[i],Pos[i],Costheta[i], Sentheta[i]);
      
    }

FuerzaElectrica(*Costheta, *Sentheta,  n_points, *Vel); 
  
  return 0;

  }

//Consideracion de la fuerza electrica para la aceleracion
int FuerzaElectrica(double *Costheta, double *Sentheta, int n_points, double *Vel){
int k;
	for(k=0; k<n_points; k++){

		if(Costheta[k]>0){
			if(Sentheta[k]<0.0001){
			Vel[k]=sqrt(Vel[k-1]*Vel[k-1]+2*d*a);
			}
		}
		else if(Costheta[k]<0){
			if(Sentheta[k]<0.0001){
			Vel[k]=sqrt(Vel[k-1]*Vel[k-1]+2*d*(-a));
			}
		}
	}
return 0;
}

//Funciones que rigen el movimiento radial

double v_prime( double v, double Vel) {
  
  return (q*v*B)/(m);
  
}

double r_prime( double v, double Vel) {
  
  return v;
  
}

//Metodo rungekutta

int rungekutta(int i, double *Vel,double *Pos,double *t, double h)
{
  
  double kx1 = v_prime(Vel[i-1],Pos[i-1]);

  double ky1 = r_prime(Vel[i-1],Pos[i-1]);
  
  
  // Paso1
  double t1 = t[i-1] + (h/2.0);
  double Vel1 = Vel[i-1] + (h/2.0) * kx1;
  double Pos1 = Pos[i-1] + (h/2.0) * ky1;
  
	
  double kx2 = v_prime(Vel1, Pos1);
  double ky2 = r_prime(Vel1, Pos1);
  
  // Paso2
  double t2 = t[i-1] + (h/2.0);
  double Vel2 = Vel[i-1] + (h/2.0) * kx2;
  double Pos2 = Pos[i-1] + (h/2.0) * ky2;
  
  
  double kx3 = v_prime(Vel2, Pos2);
  double ky3 = r_prime(Vel2, Pos2);
  
  
	// Paso3  
  double t3 = t[i-1] + h;
  double Vel3 = Vel[i-1] + h * kx3;
  double Pos3 = Pos[i-1] + h * ky3;
  
  double kx4 = v_prime(Vel3, Pos3);
  double ky4 = r_prime(Vel3, Pos3);
  
  	
  	// Paso4
  double average_kx = (1.0/6.0)*(kx1 + 2.0*kx2 + 2.0*kx3 + kx4);
  double average_ky = (1.0/6.0)*(ky1 + 2.0*ky2 + 2.0*ky3 + ky4);
  
  // Promedio de las variables
  t[i] = t[i-1] + h;
  Vel[i] = Vel[i-1] + h * average_kx;
  Pos[i] = Pos[i-1] + h * average_ky;
 
	return 0;
}


