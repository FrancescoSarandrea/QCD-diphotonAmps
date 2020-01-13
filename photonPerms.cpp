#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

//The code computes all the possible permutations of two photons which leave the order of the coloured legs unchanged. The photons are assumed to have positions 0 and 1.
using namespace std;

//std::vector<std::vector<int> > photonPerms(std::vector<int>  vec)
 int main(int argv, char ** argc) 
{

   //std::vector<std::vector<int> > matr(30); 
   int  vec[]= {0,1,2,3,4};
   int len= sizeof(vec)/sizeof(vec[0]);
   int newVec[len];


  for(int i=0; i<(len-1); ++i){
     int found1;
     
    for(int g=0; g<len; ++g){
         if(newVec[g]==0){
           found1=g;
            break;
           }
         else{
          }
         }
      swap(vec[found1], vec[i]);

    for(int m=0; m<len; ++m){
    newVec[m]=vec[m];
      }

  for(int j=0; j<(len-1); ++j){
     int found2;
    
     for(int g=0; g<len; ++g){
         if(newVec[g]==1){
           found2=g;
            break;
           }
         else{
          }
         }
     if(j==i) {
        }
     else{ swap(newVec[found2], newVec[j]);
     
    for(int k=0; k<len; ++k){
             
            if(k==0){

             cout <<"\n"<< "{" << newVec[k] << "," <<flush;
             }

             else if(k==(len-1)){

            cout << newVec[k] << "}," << endl;
             }
              else {
              
              cout  << newVec[k] << "," <<flush;
            }
       }
 // matr.push_back(newVec);
        } 
       }
   }  
  return 1;

}

 

//int main(int argc, char **argv)
//{  
   
//  int arr[6]={1,2,3,4,5,6};

 // std::vector<int> vec(arr, arr+6);

 //int  n= sizeof(vec)/sizeof(vec[0]);


  






 // std::vector<std::vector<int> > ords= photonPerms(vec);
 
 //  for(int i=0; i<30; ++i){

//      std::vector<int> row= ords[i];

//     for(int k=0; k<6; ++k){
 //            int value= row[k];
 //           if(k==0){

  //           cout <<"\n"<< value <<flush;
  //           }
  //            else{
              
    //          cout << value <<flush;
    //        }
   //    }
 //  }
//  return 1;
//}
