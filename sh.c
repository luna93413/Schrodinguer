#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>


// Definimos constantes

#define h 0.01 // Espaciado espacial
#define PI 3.141592653589793
#define N 1000 // Número de puntos espaciales
#define I _Complex_I // Definición de la unidad imaginaria
int T = 1500; // Número de pasos temporales

// Función para inicializar la función de onda Phi
void inicializar_phi(double complex *Phi, double n_ciclos) {
    double k_tilde = (2.0 * PI * n_ciclos) / N; // \tilde{k}_0 = 2π n_ciclos / N
    double x0 = (N * h) / 4.0;                 // x_0 = Nh / 4
    double sigma = (N * h) / 16.0;             // σ = Nh / 16

    for (int j = 1; j < N; j++) {
        double x = j * h; // Posición espacial
        Phi[j] = cexp(I * k_tilde * j) * cexp(-pow((x - x0), 2) / (2 * sigma * sigma)); // Φ_{j,0} = e^{i\tilde{k}_0 x} e^{-(x - x_0)^2 / (2σ^2)}
    }

    Phi[0] = 0.0;
    Phi[N] = 0.0;

    // Normalización de Phi
    double norma = 0.0;
    for (int j = 0; j <= N; j++) {
        norma += pow(cabs(Phi[j]), 2) * h; // Norma = Σ |Φ[j]|^2 * h
    }
    norma = sqrt(norma);

    for (int j = 1; j <N; j++) {
        Phi[j] /= norma; // Normalizamos Φ[j]
    }
}

// Función para inicializar el potencial V
void inicializar_potencial(double *V, double lambda, double k_tilde) {
    double V_j_h2 = lambda * k_tilde * k_tilde; // Altura del potencial
    int inicio = (2 * N) / 5;
    int fin = (3 * N) / 5;

    for (int j = 0; j <= N; j++) {
        if (j >= inicio && j <= fin) {
            V[j] = V_j_h2; // Dentro del rango del potencial
        } else {
            V[j] = 0.0; // Fuera del rango del potencial
        }
    }
}

// Función para encontrar el primer máximo local de un array PD de tamaño T
int primer_maximo_local(double *PD, int T) { 
    
    for (int n = 1; n < T - 6; n++) {
        if (PD[n] >0.01 && PD[n - 1] && PD[n] > PD[n + 5]) {
            return n;
        }
    }
    return -1; // No encontrado
}

// Calcula <j> para el caso detectada o no detectada
double valor_esperado_j(double complex *Phi, int inicio, int fin, int detectada) {
    double norma = 0.0;
    double j_esperado = 0.0;
    double complex Phi_tmp[N+1];

    // Copia Phi a Phi_tmp para no modificar el original
    for (int j = 0; j <= N; j++) {
        Phi_tmp[j] = Phi[j];
    }

    if (detectada) {
        // Anula fuera del detector
        for (int j = 0; j <= N; j++) {
            if (j < inicio || j > fin) Phi_tmp[j] = 0.0;
        }
        // Renormaliza dentro del detector
        for (int j = inicio; j <= fin; j++) {
            norma += pow(cabs(Phi_tmp[j]), 2) * h;
        }
        norma = sqrt(norma);
        if (norma == 0) return 0.0;
        for (int j = inicio; j <= fin; j++) {
            Phi_tmp[j] /= norma;
        }
        // Calcula <j>
        for (int j = inicio; j <= fin; j++) {
            j_esperado += j * pow(cabs(Phi_tmp[j]), 2) * h;
        }
    } else {
        // Anula dentro del detector
        for (int j = inicio; j <= fin; j++) {
            Phi_tmp[j] = 0.0;
        }
        // Renormaliza fuera del detector
        for (int j = 0; j <= N; j++) {
            if (j < inicio || j > fin) norma += pow(cabs(Phi_tmp[j]), 2) * h;
        }
        norma = sqrt(norma);
        if (norma == 0) return 0.0;
        for (int j = 0; j <= N; j++) {
            if (j < inicio || j > fin) Phi_tmp[j] /= norma;
        }
        // Calcula <j>
        for (int j = 0; j <= N; j++) {
            if (j < inicio || j > fin) j_esperado += j * pow(cabs(Phi_tmp[j]), 2) * h;
        }
    }
    return j_esperado;
}


int main() {
    srand((unsigned)time(NULL));

    // --- Cálculo del coeficiente de transmisión ---
    int m = 10; // Número de experimentos
    int m_T = 0;  // Número de transmisiones
    double PD[T]; // Guardar PD(n) en cada paso de tiempo
    int n_D = 0;  // Paso del primer máximo local
    double PD_nD = 0.0;

    double complex Phi[N + 1], Phi_next[N + 1], chi[N + 1];
    double V[N + 1];
    complex double a[N + 1], b[N + 1], c[N + 1], d[N + 1];
    complex double gamma[N -1]; // Coeficientes gamma

       // Parámetros iniciales
    int n_ciclos = N/16; // Número de ciclos
    double lambda = 1; // Altura del potencial
    double k_tilde = (2.0 * PI * n_ciclos) / N; // \tilde{k}_0

    // Calculamos s tilde
    double tilde_s = 1.0 / (4.0 * k_tilde * k_tilde);

    FILE *f_pd = fopen("probabilidades_pd.txt", "w");
    if (f_pd == NULL) {
    fprintf(stderr, "Error al abrir el fichero para guardar PD.\n");
    return 1;
    } 

    FILE *f_jesp = fopen("j_esperado.txt", "w");
    if (f_jesp == NULL) {
    fprintf(stderr, "Error al abrir el fichero para guardar <j>.\n");
    return 1;
    }

     FILE *f_jesp2 = fopen("j_esperado2.txt", "w");
    if (f_jesp2 == NULL) {
    fprintf(stderr, "Error al abrir el fichero para guardar <j>.\n");
    return 1;
    }

    for (int experimento = 0; experimento <=m; experimento++) {
    // Declaramos las variables

    a[0] = 0.0; // a_0 = 0
    b[0] = 0.0; // b_0 = 0
    c[0] = 0.0; // c_0 = 0
    d[0] = 0.0; // d_0 = 0

    // Inicializamos la función de onda y el potencial
    inicializar_phi(Phi, n_ciclos);
    inicializar_potencial(V, lambda, k_tilde);

     // Abrimos el fichero para guardar las normas
     FILE *f_normas = fopen("normas.txt", "w");
     if (f_normas == NULL) {
         fprintf(stderr, "Error al abrir el fichero para guardar las normas.\n");
         return 1;
     }

     // Abrimos el fichero para guardar los datos de la función de onda
    FILE *f_datos = fopen("schrodinger.dat", "w");
    if (f_datos == NULL) {
        fprintf(stderr, "Error al abrir el fichero para guardar los datos.\n");
        return 1;
    }


      for (int j = 1; j < N; j++) {
        a[j] = 1.0;
        b[j] = -2.0 + 2.0 * I / tilde_s - V[j];
        c[j] = 1.0;
      }

      // Calculamos alpha fuera del bucle temporal
      complex double alpha[N];
      alpha[N-1] = 0.0;

        for (int j = N-1; j > 0; j--) {
            gamma[j] = b[j] + c[j] * alpha[j];
            alpha[j-1] = -a[j] / gamma[j];
        }

        int j;
    // Iteración en el tiempo
    
    for (int n = 0; n < T; n++) {

         // Guardar datos para animacion.py: x, |Phi|^2, V
        for (j = 0; j <= N; j++) {
            double x = j * h;
            double densidad = pow(cabs(Phi[j]), 2);
            fprintf(f_datos, "%.8f,%.8f,%.8f\n", x, densidad, V[j]);
        }
        fprintf(f_datos, "\n"); // Línea en blanco entre bloques

        // Construimos los coeficientes del sistema tridiagonal
        for (int j = 1; j < N; j++) {
            d[j] = 4.0 * I / tilde_s * Phi[j];
        }

        // Resolvemos el sistema tridiagonal para chi
        complex double beta[N]; // Beta aún depende del tiempo
        beta[N - 1] = 0.0 + 0.0*I;      // β_{N-1,n} = 0
        // Orden descendente
        for (int j = N - 1; j > 0; j--) {
            beta[j - 1] = (d[j] - c[j] * beta[j]) / gamma[j];
        }

        // Orden ascendente para calcular χ
        chi[0] = 0.0;
        chi[N] = 0.0;
        for (int j = 0; j < N - 1; j++) {
            chi[j + 1] = alpha[j]*chi[j] + beta[j];
        }

        // Calculamos Phi_next
        for (int j = 1; j < N; j++) {
            Phi_next[j] = chi[j] - Phi[j];
        }

        // Actualizamos Phi
        for (int j = 1; j < N; j++) {
            Phi[j] = Phi_next[j];
        }

        // Calculamos la norma de la función de onda
        double norma = 0.0;
        for (int j = 1; j <N; j++) {
            norma += pow(cabs(Phi[j]), 2) * h; // Norma = Σ |Φ[j]|^2 * h
        }

        // Guardamos la norma en el fichero
        fprintf(f_normas, "Paso %d: Norma = %.10f\n", n, norma);
        

        //COEFICIENTE DE TRANSMISIÓN
        if(experimento==0) {
            double pd = 0.0;
            for (int j = (4*N)/5; j <= N; j++) {
                pd += pow(cabs(Phi[j]), 2) * h;
            }
            PD[n] = pd;

        
             fprintf(f_pd, "%d %.20f\n", n, PD[n]);

            double j_esperado = 0.0;
            for (int j = 0; j <= N; j++) {
                 j_esperado += j * pow(cabs(Phi[j]), 2) * h;
             }
             fprintf(f_jesp, "%d %.10f\n", n, j_esperado);
        }
   
    }
        fclose(f_normas);
        fclose(f_datos);

        // Buscar el primer máximo local de PD(n)
        if(experimento==0) {
            n_D = primer_maximo_local(PD, T);
          printf("Primer máximo local de PD(n) en n = %d\n", n_D);

           PD_nD = PD[n_D];
             printf("PD(n_D) = %.20f\n", PD_nD);
        }
       
        T=n_D;
        int detectada;
       if(experimento>0) {
        double p = (double)rand() / ((double)RAND_MAX + 1.0);
        printf("Número aleatorio p = %.20f\n", p);
             if (p < PD_nD) {
                 m_T++;
                 detectada=1;
            }else 
            {
                 detectada=0;
             }
           fprintf(f_jesp2, "%d %.10f %d\n", experimento, valor_esperado_j(Phi, (4*N)/5, N, detectada), m_T);
     }
    

}
     fclose(f_pd);
     fclose(f_jesp);
     fclose(f_jesp2);
    double K = (double)m_T / m;
    printf("Número de transmisiones m_T = %d\n", m_T);
    printf("Coeficiente de transmisión K = %.6f\n", K);
    // --- Fin del cálculo de transmisión ---

    return 0;
}


   
    
        