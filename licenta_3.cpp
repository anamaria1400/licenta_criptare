#include <iostream>
#include <cmath>
#include <math.h>
#include <string.h>
#include <list>
#include <iterator>
using namespace std;

class CommonMethods {
public:
    int gen = 7; //initializarea generatorului
    double SCALE_FACTOR = 100000; //initializarea factorului de scala

    //FUNCTIA SIN
    double sinFunction(double k, double l) {
        return sin((k - l) / 2);
    }
    //GASIREA CUTIEI
    int findBoxNumber(double sinValue) {
        if (sinValue >= -1.0 && sinValue < -0.8) { return 1; }
        else if (sinValue >= -0.8 && sinValue < -0.6) { return 2; }
        else if (sinValue >= -0.6 && sinValue < -0.4) { return 3; }
        else if (sinValue >= -0.4 && sinValue < -0.2) { return 4; }
        else if (sinValue >= -0.2 && sinValue < 0.0) { return 5; }
        else if (sinValue >= 0.0 && sinValue < 0.2) { return 6; }
        else if (sinValue >= 0.2 && sinValue < 0.4) { return 7; }
        else if (sinValue >= 0.4 && sinValue < 0.6) { return 8; }
        else if (sinValue >= 0.6 && sinValue < 0.8) { return 9; }
        else if (sinValue >= 0.8 && sinValue < 1.0) { return 10; }
        else return 0;
    }
    //BITUL CORESPUNZATOR NUMARULUI CUTIEI
    int findByteNumber(int boxNumber) {
        return (int)fmod((pow(gen, boxNumber)), 255);
    }
    //XOR INTRE 2 NUMERE
    int* bitXor(int* k, int* l) {
        int* vect_xor = new int[8]();
        for (int i = 0; i < 8; i++) {
            vect_xor[i] = k[i] ^ l[i];
        }
        return vect_xor;
    }
    //TRANSFORMA NUMARUL LONG IN VECTOR
    int* convertVect(long n) {
        int* vect = new int[8]();
        for (int i = 7; i >= 0; i--) {
            vect[i] = n % 10;
            n /= 10;
        }
        return vect;
    }
    //CONVERTESTE NR DIN ZECIMAL IN BINAR
    long convertBinary(int n) {
        int remainder;
        long binary = 0, i = 1;

        while (n != 0) {
            remainder = n % 2;
            n = n / 2;
            binary = binary + (remainder * i);
            i = i * 10;
        }

        return binary;
    }
    //CONVERTESTE DIN BINAR IN ZECIMAL
    int binaryToDecimal(int n) {
        int num = n;
        int dec_value = 0;
        int base = 1;// Initializeaza valoarea bazei cu 1, adica 2^0
        int temp = num;
        while (temp) {
            int last_digit = temp % 10;
            temp = temp / 10;

            dec_value += last_digit * base;

            base = base * 2;
        }
        return dec_value;
    }

    //ALEGE PERMUTAREA IN FUNCTIE DE NUMARUL CUTIEI
    /*int** findPermutationNumber(int boxNumber) {
        list<int[8][8]> permutations; //creeaza un ArrayList         
        int matrix1[8][8] = {{0,0,0,0,0,0,0,1}, {0,0,0,0,0,0,1,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,1,0,0,0,0}, {0,0,1,0,0,0,0,0}, {0,1,0,0,0,0,0,0}, {1,0,0,0,0,0,0,0}};
        int matrix2[8][8] = {{1,0,0,0,0,0,0,0}, {0,1,0,0,0,0,0,0}, {0,0,1,0,0,0,0,0}, {0,0,0,1,0,0,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,0,0,1,0}, {0,0,0,0,0,0,0,1}};
        int matrix3[8][8] = {{0,1,0,0,0,0,0,0}, {1,0,0,0,0,0,0,0}, {0,0,1,0,0,0,0,0}, {0,0,0,1,0,0,0,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,0,0,0,1,0}, {0,0,0,0,0,0,0,1}};
        int matrix4[8][8] = {{0,1,0,0,0,0,0,0}, {0,0,1,0,0,0,0,0}, {0,0,0,1,0,0,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,0,0,1,0}, {0,0,0,0,0,0,0,1}, {1,0,0,0,0,0,0,0}};
        int matrix5[8][8] = {{0,0,0,0,0,0,1,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,1,0,0,0,0}, {0,0,1,0,0,0,0,0}, {0,1,0,0,0,0,0,0}, {1,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,1}};
        int matrix6[8][8] = {{0,0,0,0,0,0,0,1}, {0,0,0,0,0,1,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,0,0,0,1,0}, {0,0,0,1,0,0,0,0}, {0,1,0,0,0,0,0,0}, {0,0,1,0,0,0,0,0}, {1,0,0,0,0,0,0,0}};
        int matrix7[8][8] = {{0,0,1,0,0,0,0,0}, {0,1,0,0,0,0,0,0}, {1,0,0,0,0,0,0,0}, {0,0,0,1,0,0,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,0,0,0,1}, {0,0,0,0,0,0,1,0}};
        int matrix8[8][8] = {{0,0,0,0,0,0,0,1}, {0,0,0,0,0,0,1,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,1,0,0,0}, {0,0,1,0,0,0,0,0}, {0,1,0,0,0,0,0,0}, {1,0,0,0,0,0,0,0}, {0,0,0,1,0,0,0,0}};
        int matrix9[8][8] = {{0,0,0,0,0,1,0,0}, {0,0,0,0,0,0,0,1}, {0,0,0,0,0,0,1,0}, {0,0,0,0,1,0,0,0}, {0,0,0,1,0,0,0,0}, {0,0,1,0,0,0,0,0}, {0,1,0,0,0,0,0,0}, {1,0,0,0,0,0,0,0}};
        int matrix10[8][8] = {{0,0,0,0,0,0,0,1}, {0,0,0,0,0,0,1,0}, {0,0,0,0,0,1,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,1,0,0,0,0}, {1,0,0,0,0,0,0,0}, {0,0,1,0,0,0,0,0}, {0,1,0,0,0,0,0,0}};
        permutations.push_back(matrix1); permutations.push_back(matrix2);
        permutations.push_back(matrix3); permutations.push_back(matrix4);
        permutations.push_back(matrix5); permutations.push_back(matrix6);
        permutations.push_back(matrix7); permutations.push_back(matrix8);
        permutations.push_back(matrix9); permutations.push_back(matrix10); //adauga matricile in "permutations"     
        int a = boxNumber - 1;
        list<int[8][8]>:: iterator it = permutations.begin();
        advance(it, a);
        
        return (int**)*it; //permutarea finala ia valoarea permutarii corespunzatoare     
    }*/
    //INMULTIRE INTRE DOUA MATRICI
    /*double** multiplyMatrices(double m[3][3], double n[3][3]) {
        double** product = new double*;
        for (int i = 0; i < 8; i++) { //se parcurg liniile primei matrici [r_m]
            for (int k = 0; k < 8; k++) { //se parcurg coloanele primei matrici [c_m]
                for (int j = 0; j < 8; j++) { //se parcurg coloanele celei de-a doua matrici [c_n]
                    product[i][j] += m[i][k] * n[k][j];
                }

            }
        }
        return product;
    }*/
    //INMULTIRE INTRE O MATRICE SI UN VECTOR
    /*int* multiply(int** m, int* n) {
        int** a = new int* [8];
        for (int i = 0; i < 8; ++i)
            a[i] = new int[8];
        a = m;
        int* product = new int [8]();
        for (int i = 0; i < 8; i++) // i exceeds the index
        {
            product[i] = 0; // warning C6200
            
        }
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                product[j] += (a[j][k] * n[k]);
            }
        }
        return product;
    }
    */
    //TRANSFORMA VECTORUL IN NUMAR
    int convertVectInNumber(int* vector) {
        int* vect = vector;
        int number = 0;
        for (int i = 0; i < 8; i++) {
            number = number*10 + vect[i];
        }
        return (int)number;
    }
};

class Encryption : public CommonMethods {
public:
    double x0, x1, x2, a, b;
    double m;

    Encryption(double x0_0, double x1_0, double x2_0, double a_0, double b_0) {
        x0 = x0_0; x1 = x1_0; x2 = x2_0; a = a_0; b = b_0;
    }

    //ADAUGA CIFRUL SI TEXTUL CRIPTAT IN CATE UN VECTOR SI INTOARCE CIFRUL
    list<double> encrypt(list<int> input, list<int> cipherText)//rezolva sistemul Henon3D si gaseste cifrul
    {
        list<double> cipher; //initializeaza un ArrayList
        cipher.push_back(x2); //prima cheie este al treilea termen din cond initiale
        for (int i = 0; i <= input.size(); i++) {
            m = resultAfterEncryption(i);
            cipherText.push_back((int)m); //adauga in vectorul cipherText textul criptat
            double x0_new = a - x1 * x1 - b * x2 + m / SCALE_FACTOR;
            double x1_new = x0;
            double x2_new = x1; //sistemul Henon3D
            x0 = x0_new;
            x1 = x1_new;
            x2 = x2_new; //se actualizeaza valorile calculate cu sistemul Henon 3D
            cipher.push_back(x2);
        }
        cipher.push_back(x1);
        cipher.push_back(x0);
        return cipher;

    }
    //INTOARCE VALOARE IN ZECIMAL
    int resultAfterEncryption(int bt) //criptarea textului in clar prin alg propus
    {
        double sinValue = sinFunction(x1, x2);
        int boxNumber = findBoxNumber(sinValue);
        int inputValue = bt;
        int byteNumber = findByteNumber(boxNumber);
        //int** permutationNumber = findPermutationNumber(boxNumber);
        int permutation[8][8] = { {0,0,0,0,0,0,0,1}, {0,0,0,0,0,1,0,0}, {0,0,0,0,1,0,0,0}, {0,0,0,0,0,0,1,0}, {0,0,0,1,0,0,0,0}, {0,1,0,0,0,0,0,0}, {0,0,1,0,0,0,0,0}, {1,0,0,0,0,0,0,0} };
        long convert_bt = convertBinary(bt);
        long convert_byte = convertBinary(byteNumber);
        int* binaryBitXor = bitXor(convertVect(convert_bt), convertVect(convert_byte));
        //int* resultMultiplyMatrices = multiply((int**)permutation, binaryBitXor);
        int* product = new int[8]();

        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                product[j] += (permutation[j][k] * binaryBitXor[k]);
            }
        }

        int x = convertVectInNumber(product);
        return binaryToDecimal(x);
    }
    //CALCULEAZA EXPONENTII
    void LyapunovExponents()//calculeaza exponentul Lyapunov
    {
        const int N = 3; //ordinul sistemului Henon
        int iter = 5000; //numarul de iteratii
        double lambda[N] = { 0 }; //initializeaza vectorul in care se salveaza exponentii
        double norm[N] = { 0 }; //vectorul in care se vor pune normarile
        double v[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} }; //se defineste v ca matrice unitate
        double J[3][3] = { {0,-2 * x1, -b},{1,0,0},{0,1,0} }; //se defineste Jacobianului
        double w[3][3] = { {0} };
        for (int i = 0; i < iter; i++) {
            //double** w = multiplyMatrices(v, J); //se face inmultirea dintre Jacobian si matricea unitate, formula (2.2)
           
            for (int d = 0; d < N; d++) { //se parcurg liniile primei matrici [r_m]
                for (int e = 0; e < N; e++) { //se parcurg coloanele primei matrici [c_m]
                    for (int f = 0; f < N; f++) { //se parcurg coloanele celei de-a doua matrici [c_n]
                        w[d][f] += v[d][e] * J[e][f];
                    }
                }
            }
            // ortonormare prin metoda Gram-Schmidt pentru prima linie a lui v
            norm[0] = sqrt(w[0][0] * w[0][0] + w[0][1] * w[0][1] + w[0][2] * w[0][2]); //prima normare pentru prima linie a lui v
            for (int j = 0; j < N; j++) {
                 v[0][j] = w[0][j] / norm[0]; //prima linie a lui v se normeaza la prima valoare a lui norm
            }

            // ortonormare prin metoda Gram-Schmidt pentru a doua linie a lui v
            for (int j = 0; j < N; j++) {
                v[1][j] = w[1][j] - (w[1][0] * v[0][0] + w[1][1] * v[0][1] + w[1][2] * v[0][2]) * v[0][j]; //scrierea lui v ca numaratorul fractiei (formula (2.4))
            }
            norm[1] = sqrt(v[1][0] * v[1][0] + v[1][1] * v[1][1] + v[1][2] * v[1][2]); //a doua normare pentru a doua linie a lui v
            for (int j = 0; j < N; j++) {
                v[1][j] = v[1][j] / norm[1]; //a doua linie a lui v se normeaza la a doua valoare a lui norm
            }

            // ortonormare prin metoda Gram-Schmidt pentru a treia linie a lui v
            for (int j = 0; j < N; j++) {
                v[2][j] = w[2][j] - (w[2][0] * v[1][0] + w[2][1] * v[1][1] + w[2][2] * v[1][2]) * v[1][j] - (w[2][0] * v[0][0] + w[2][1] * v[0][1] + w[2][2] * v[0][2]) * v[0][j];  //scrierea lui v ca numaratorul fractiei (formula (2.4))
            }
            norm[2] = sqrt(v[2][0] * v[2][0] + v[2][1] * v[2][1] + v[2][2] * v[2][2]); //a treia normare pentru a treia linie a lui v
            for (int j = 0; j < N; j++) {
                v[2][j] = v[2][j] / norm[2]; //a treia linie a lui v se normeaza la a treia valoare a lui norm
            }              
            for (int j = 0; j < N; j++) { //iterare dupa numarul N
                lambda[j] += log(norm[j]); //se aduna la fiecare lambda valoarea normata conform formulei (2.1)
            }
            double x0_new = a - x1 * x1 - b * x2;
            double x1_new = x0;
            double x2_new = x1; //sistemul Henon 3D
            x0 = x0_new;
            x1 = x1_new;
            x2 = x2_new; //se actualizeaza valorile calculate cu sistemul Henon 3D
            J[0][1] = (-2) * x1; //se actualizeaza termentul Jacobianului care il contine pe x1 cu noua valoare
            
            for (int j = 0; j < N; j++) {
                lambda[j] = lambda[j] / (log(2) / log(2)) / iter; //formula (2.1) -> log2(2)
            }
            cout << lambda[0] << " " << lambda[1] << " " << lambda[2]; //se scriu in fisier exponentii Lyapunov

        }
    }

    //CALCULEAZA VALORILE PENTRU DIAGRAMA DE BIFURCATIE
    void bifurcationDiagram()//calculeaza valorile pentru diagrama de bifurcatie
    {
        const int iter = 5000; //numarul de iteratii
        double a; //se declara parametrul de bifurcatie a
        pair<double, list<double>> bifurcation; //se creeaza o mapa pentru a corela o valoare a lui a cu toate valorile lui x0 pentru acel a
        for (a = 0.001; a < 2; a += 0.001)
        {
            double x0_l[iter + 1];
            double x1_l[iter + 1];
            double x2_l[iter + 1]; //se creeaza trei liste pentru x0, x1 si x2
            x0_l[0] = x0;
            x1_l[0] = x1;
            x2_l[0] = x2; //se adauga in liste x0, x1 si x2
            for (int j = 1; j < iter; j++) {
                x0_l[j] = a - x1_l[j - 1] * x1_l[j - 1] - b * x2_l[j - 1];
                x1_l[j] = x0_l [j - 1] ;
                x2_l[j] = x1_l[j - 1]; //se calculeaza x0, x1 si x2 cu sistemul Henon 3D
            }
            /*for (int j = 0; j < iter - 2000; j++)
            {
                x0_l.remove(j); //se sterg primele 2000 de valori ale lui x0 (pana se stabilizeaza sistemul)
                j--; //se scade j pentru ca, daca am sters un element, urmatorul se incarca pe pozitia curenta si la urmatorul j se sare peste el
            }*/
            //bifurcation.insert({ a, x0_l }); //se adauga perechea a si vectorul x0 in mapa
            for (int j = 2000; j <= 5000; j++) {
                cout << x0_l[j] << " ";
            }
        }
        /*for (Entry<Double, ArrayList<Double>> entry : bifurcation.entrySet()) { //iterare dupa mapa
            for (int j = 0; j < entry.getValue().size(); j++) { //iterare dupa valorile lui x0 asociate unui a
                writeBifurcation.write(entry.getValue().get(j) + " "); //se scriu valorile in fisier
            }
            cout << endl;
        }*/
       
    }
    
};


int main() {

    double x0, x1, x2, a, b;
    cout << "Introduceti cheia:" << endl;
    cout << "x0 = "; cin >> x0; cout << endl;
    cout << "x1 = "; cin >> x1; cout << endl;
    cout << "x2 = "; cin >> x2; cout << endl;
    cout << "a = "; cin >> a; cout << endl;
    cout << "b = "; cin >> b; cout << endl;
    
    Encryption enc1(x0, x1, x2, a, b);

    cout << enc1.resultAfterEncryption((int)'E') << endl;
    cout << enc1.resultAfterEncryption((int)'t') << endl;
    cout << enc1.resultAfterEncryption((int)'t') << endl;
    cout << enc1.resultAfterEncryption((int)'I') << endl;


    list<double> cipher; //se creeaza vectorul in care se adauga valorile lui x2 (cifrul)         
    list<int> cipherText; //se creeaza vectorul in care se adauga textul criptat
    list<int> input = { 'E', 't', 't','I' };
    for (double d : enc1.encrypt(input, cipherText)) { //iterare dupa  valorile cifrului                          
        cipher.push_back(d); //se adauga valorile cifrului in vectorul dedicat                       
        cout << d << endl;
    }        
    
    enc1.LyapunovExponents();
    enc1.bifurcationDiagram();
    return 0;
}
