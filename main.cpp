#include <iostream>
#include <Eigen/Dense>
#include <iomanip> 
#include <math.h>
#include <fstream>
#include <vector>

using namespace Eigen;

struct datapoint{
    float x;
    float y;
};

float myround(float);
void processmoddatafile(std::ifstream& , std::vector<float>& ); 
void processorderedmoddatafile(std::ifstream&, std::vector<float>&);
float matrixcondition(MatrixXd);
std::vector<float> divdiff(std::vector<float>);
float positive(float);
void InsertionSort(std::vector<datapoint>&);

int main(int argc, char** argv)
{   
    std::ifstream mod1datafile, mod2datafile, mod3datafile;
    std::vector<float> mod1datafilecont, mod2datafilecont;
	if(argc != 6){
		std::cout << "You need to supply 5 arguments to this program." << "\n";
		std::cout << "The number of arguments you supplied was " << argc - 1 << ".\n";
		return -1;
	}
    mod1datafile.open(argv[1]);
    if (!mod1datafile.is_open())
      std::cout << "Could not open file\n";
    else {
        processmoddatafile(mod1datafile, mod1datafilecont);
    }
    int mod1funcdeg = atoi(argv[2]);
    MatrixXd a(mod1funcdeg+1, mod1datafilecont.size()/2);
    for(int mod1data2mat0 = 0; mod1data2mat0 <=  mod1funcdeg;mod1data2mat0++ ){
        for(int mod1data2mat1 = 0; mod1data2mat1 < (mod1datafilecont.size()/2); mod1data2mat1++){
            a(mod1data2mat0,mod1data2mat1) = pow(mod1datafilecont[mod1data2mat1*2],mod1data2mat0);
        }
    }
    MatrixXd b(1, mod1datafilecont.size()/2);
    for(int mod1data2mat2 = 0; mod1data2mat2 < mod1datafilecont.size()/2;mod1data2mat2++){
        b(0,mod1data2mat2) = mod1datafilecont[2 * mod1data2mat2 + 1];
    }
    MatrixXd aa(mod1funcdeg+1,mod1funcdeg+1);
    MatrixXd bb(1,mod1datafilecont.size()/2);
    VectorXd c(mod1funcdeg+1);
    VectorXd d(mod1funcdeg+1);
    aa = a * a.transpose();       //A
    bb = b * a.transpose();       //B in matrix form
    //std::cout << aa << "\n";  //debug
    //std::cout << bb << "\n";  //debug
    for(int mod1data2mat3 = 0; mod1data2mat3 < mod1funcdeg+1; mod1data2mat3++){
        c[mod1data2mat3] = bb(0,mod1data2mat3);
    }
    d = aa.inverse() * c;         //C
    for(int roundcount = 0;roundcount < mod1funcdeg+1;roundcount++){
        d(roundcount) = myround(d(roundcount));
    }
    std::cout << "Co-Efficients for mod 1 (x^0,x^1,...,x^n):"<< "\n" <<  d << "\n";
    std::cout << "Condition for Matrix A:" << matrixcondition(aa) << "\n";
       
    mod2datafile.open(argv[3]);
    if (!mod2datafile.is_open())
      std::cout << "Could not open file\n";
    else {
        processmoddatafile(mod2datafile, mod2datafilecont);
    }
    int mod2funcdeg = atoi(argv[4]);
    MatrixXd mod2prea(mod2funcdeg+1, mod2datafilecont.size()/2);
    for(int mod2data2mat0 = 0; mod2data2mat0 <=  mod2funcdeg;mod2data2mat0++ ){
        for(int mod2data2mat1 = 0; mod2data2mat1 < (mod2datafilecont.size()/2); mod2data2mat1++){
            if (mod2data2mat0 == 0){
                mod2prea(mod2data2mat0,mod2data2mat1) = 1;
            }
            else if (mod2data2mat0 == 1){
                mod2prea(mod2data2mat0,mod2data2mat1) = mod2datafilecont[mod2data2mat1*2];
            }
            else{
                mod2prea(mod2data2mat0,mod2data2mat1) = ((2*mod2data2mat0-1)*mod2prea(1,mod2data2mat1)*mod2prea(mod2data2mat0-1,mod2data2mat1)-((mod2data2mat0 - 1)*mod2prea(mod2data2mat0-2,mod2data2mat1)))/(mod2data2mat0);
            }
        }
    }
    //std::cout << mod2prea << "\n";    //debug
    MatrixXd mod2preb(1, mod2datafilecont.size()/2);
    for(int mod2data2mat2 = 0; mod2data2mat2 < mod2datafilecont.size()/2;mod2data2mat2++){
        mod2preb(0,mod2data2mat2) = mod2datafilecont[2 * mod2data2mat2 + 1];
    }
    //std::cout << mod2preb << "\n";    //debug
    MatrixXd mod2a(mod2funcdeg+1,mod2funcdeg+1);
    MatrixXd mod2b(1,mod2datafilecont.size()/2);
    VectorXd mod2c(mod2funcdeg+1);
    VectorXd mod2d(mod2funcdeg+1);
    mod2a = mod2prea * mod2prea.transpose();       //A
    mod2b = mod2preb * mod2prea.transpose();       //B in matrix form
    //std::cout << mod2a << "\n";  //debug
    //std::cout << mod2b << "\n";  //debug
	//std::cout << mod2a.inverse() << "\n";		//debug
    for(int mod2data2mat3 = 0; mod2data2mat3 < mod2funcdeg+1; mod2data2mat3++){
        mod2c[mod2data2mat3] = mod2b(0,mod2data2mat3);
    }    
    mod2d = mod2a.inverse() * mod2c;         //C
    for(int roundcount1 = 0;roundcount1 < mod2funcdeg+1;roundcount1++){
        mod2d(roundcount1) = myround(mod2d(roundcount1));
    }
    std::cout << "Co-Efficients for mod 2 (Legendre Polynomials):"<< "\n" <<  mod2d << "\n";
    std::cout << "Condition for Matrix A:" << matrixcondition(mod2a) << "\n";
    
     std::vector<float> mod3datacont, answer;
    mod3datafile.open(argv[5]);
    if (!mod3datafile.is_open())
      std::cout << "Could not open file\n";
    else {
        processorderedmoddatafile(mod3datafile, mod3datacont);
    }   
    answer = divdiff(mod3datacont);
    std::cout << "Co-Efficients for mod 3 (divided-difference):" << "\n"; 
    for (int count0 = 0; count0 < answer.size(); count0++) {
        std::cout << myround(answer[count0]) << '\n';
    }
    
    return 0;
}

float myround(float a){
    a = a*10000+.5;
    a = floor(a);
    a = a/10000;
    return a;                 
}

void processmoddatafile(std::ifstream& a, std::vector<float>& b){
    float c;
    while(a >> c){
        b.push_back(c);
    }
}

void processorderedmoddatafile(std::ifstream& a, std::vector<float>& b){
    float c;
    datapoint data;
    std::vector<datapoint> d;
    while(a >> c){
        b.push_back(c);
    }
    for(int ordered0 = 0; ordered0 < b.size(); ordered0++){
        if(ordered0 % 2 == 0){
            data.x = b[ordered0];
        }
        if(ordered0 % 2 == 1){
            data.y = b[ordered0];
            d.push_back(data);
            data = datapoint();            
        }
    }
    InsertionSort(d);
    b.clear();
    for(int ordered1 = d.size()-1; ordered1 > -1; ordered1--){
        b.push_back(d[ordered1].x);
        b.push_back(d[ordered1].y);
    }
}

float matrixcondition(MatrixXd matrix){
    float a=0,b=0,c=0,d=0;
    MatrixXd inv(matrix.rows(),matrix.cols());
    inv = matrix.inverse();
    for(int matrixcondition0 = 0; matrixcondition0 < matrix.rows(); matrixcondition0++){
        for(int matrixcondition1 = 0; matrixcondition1 < matrix.cols(); matrixcondition1++){
            c = c + positive(matrix(matrixcondition0,matrixcondition1));
            d = d + positive(inv(matrixcondition0,matrixcondition1));
        }
        if(c>a){
            a = c;
        }if(d>b){
            b = d;
        }
        c = 0;
        d = 0;
    }
a = a * b;
return a;
}

std::vector<float> divdiff(std::vector<float> a) {
    std::vector<float> answer;
    int datanum = a.size();
    int newdata = 0;
    float maxadd = 0;
    for (int newdataconst = 0; newdataconst < a.size() / 2; newdataconst++) {
        newdata = newdata + newdataconst;
        maxadd = float (newdataconst);
    }
    answer.push_back(a[1]);
    for (int stage1 = 0; stage1 < (datanum / 2) - 1; stage1++) {
        a.push_back((a[2 * stage1 + 3] - a[2 * stage1 + 1]) / (a[2 * stage1 + 2] - a[2 * stage1]));
        if(stage1 == 0){
           answer.push_back((a[2 * stage1 + 3] - a[2 * stage1 + 1]) / (a[2 * stage1 + 2] - a[2 * stage1])); 
        }
    }
    float iterator = 2;
    float counter = 0;
    int counteriterator = (datanum / 2) - 2;
    int sizecounter = a.size();
    if (datanum / 2 > 2) {
        for (int stage2 = datanum + 1; stage2 < datanum + newdata - 1; stage2++) {
            a.push_back((a[stage2] - a[stage2 - 1]) / (a[2 * counter + 2 * iterator] - a[2 * counter + 0]));      
            if (counter == 0){				
                answer.push_back((a[stage2] - a[stage2 - 1]) / (a[2 * counter + 2 * iterator] - a[2 * counter + 0]));
            }
            counter++;
            if (a.size() == sizecounter + counteriterator) {
				sizecounter = a.size();
                counter = 0;
                counteriterator--;
                iterator++;
                stage2++;
            }

        }
    }
    return answer;
}

float positive(float number){
	if(number < 0){
		number = number * -1;
	}
	return number;
}

void InsertionSort( std::vector<datapoint> &num)
{
     int i, j, numLength = num.size();
     datapoint key;
     for(j = 1; j < numLength; j++)    
    {
           key = num[j];
           for(i = j - 1; (i >= 0) && (num[i].x < key.x); i--)   
          {
                 num[i+1] = num[i];
          }
         num[i+1] = key;    
     }
     return;
}




