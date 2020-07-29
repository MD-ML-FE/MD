//aauthor: llw 2020_7_29
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <math.h>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace Eigen;
using namespace std;

void MIN_MAX_CENTRE ( double *distance_BF, double *distance_maxF, double *distance_minF , int *d_min_index , int *d_max_index , int coutttt_B , double *centre_x , double *centre_y , double *centre_z , int crystal_size_index[] , double x[] ,double y[] , double z[] ) {
	
	*distance_BF = 0 ;
  *distance_maxF = 0;
  *distance_minF = 10000;
  *d_min_index = 100000;
  *d_max_index = 100000;
  for ( int i_max = 0 ; i_max < coutttt_B ; i_max ++ ){
		*distance_BF = ( ( *centre_x - x[i_max] ) * ( *centre_x - x[i_max] ) ) + ( ( *centre_y - y[i_max] ) * ( *centre_y - y[i_max] ) ) + ( ( *centre_z - z[i_max] ) * ( *centre_z - z[i_max] ) ) ;
		if ( *distance_minF > *distance_BF ) {
	  	*distance_minF = *distance_BF ;
	  	*d_min_index = i_max;
		}
		if ( *distance_maxF < *distance_BF ) {
	  	*distance_maxF = *distance_BF ;
	  	*d_max_index = i_max;
		}
  } 
  cout<<"min: "<<*distance_minF<<" "<<sqrt(*distance_minF)<<" ["<<*d_min_index<<"] "<<"max: "<<*distance_maxF<<" "<<sqrt(*distance_maxF)<<" ["<<*d_max_index<<"] "<<endl;
	
}

void CENTRE ( int coutttt_B , double *centre_x , double *centre_y , double *centre_z , double x[] ,double y[] , double z[] ) {
	*centre_x = 0 ;
	*centre_y = 0 ;
	*centre_z = 0 ;
	for ( int i_max = 0 ; i_max < coutttt_B ; i_max ++ ) {
		*centre_x = *centre_x + x[i_max] ;
		*centre_y = *centre_y + y[i_max] ;
		*centre_z = *centre_z + z[i_max] ;
	}
	*centre_x = *centre_x / coutttt_B ;
	*centre_y = *centre_y / coutttt_B ;
	*centre_z = *centre_z / coutttt_B ;
}

int main()
{
    cout << "Hello world!" << endl;
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + (0,1);
    std::cout << m << std::endl;
    //
    Eigen::Vector3d v_3d;
    v_3d << 1,2,3;
    Eigen::Matrix3d matrix_33 = Eigen::Matrix3d::Zero();
    //
    cout<<"matrix_33:\n"<<matrix_33<<endl;
    cout<<"v_3d:\n"<<v_3d<<endl;
    cout<<endl;
    cout<<"zhuanzhi:\n"<<v_3d.transpose()<<endl;
    cout<<"neijiA:\n"<<v_3d * (v_3d.transpose())<<endl;
    matrix_33=v_3d * (v_3d.transpose());
    matrix_33 << 4,0,0,0,3,1,0,1,3;
    cout<<"neijiB:\n"<<2*(v_3d.transpose())*v_3d<<endl;
    cout<<"matrix_33 :\n"<<matrix_33<<endl;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>eigen_solver_cs(matrix_33.transpose()*matrix_33);
    cout << "Eigen value = \n" << eigen_solver_cs.eigenvalues() << endl;
    cout << "Eigen vectors = \n" << eigen_solver_cs.eigenvectors() << endl;
    //
    Matrix3d A;
    A << 1, -1, 1, 0, 2, -3, 0, 0, 1;
    cout << "Here is a 3x3 matrix, A:" << endl << A << endl << endl;
    EigenSolver<Matrix3d> es(A);
    Matrix3d D = es.pseudoEigenvalueMatrix();
    Matrix3d V = es.pseudoEigenvectors();
    cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
    //read_gro
    cout.precision(16);
    ifstream OpenFile("panding.gro", ios::in);
    if (OpenFile.fail())
    {
	 cout<<"Can not open target file"<<endl;
         return 0;
    }
    std::string lineStr;
    std::vector<int> Xu ;
    std::vector<string> typeXu ;
    std::vector<double> x ;
    std::vector<double> y ;
    std::vector<double> z ;
    double boxx;
    double boxy;
    double boxz;
    int ii = 0;
    int b_N = 0 ;
    if (OpenFile)
    {
	int i = 0;
	while (getline(OpenFile,lineStr))
	{
	    i++;
	    if (i==2) {
		//
		string bN_0;
		istringstream bs(lineStr);
		bs>>bN_0;
		b_N = atoi(bN_0.c_str());
		cout << b_N <<endl;
	    }
	    if ( (i >= 3) && (i < (3+b_N)) ) {
		ii ++;
		string a1,a2,a3,a4,a5,a6;
		a1 = lineStr.substr(0,8);
		a2 = lineStr.substr(9,6);//type
		a3 = lineStr.substr(15,5);//XU
		a4 = lineStr.substr(21,8);
		a5 = lineStr.substr(29,8);
		a6 = lineStr.substr(37,8);
		//cout<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a5<<","<<a6<<endl;
		Xu.push_back(atoi(a3.c_str()));
		typeXu.push_back(a2.c_str());
		x.push_back(atof(a4.c_str()));
		y.push_back(atof(a5.c_str()));
		z.push_back(atof(a6.c_str()));
	    }
	    if (i == (3+b_N)) {
		string box_x,box_y,box_z;
		istringstream ip(lineStr);
		ip>>box_x>>box_y>>box_z;
		boxx = atof(box_x.c_str());
		boxy = atof(box_y.c_str());
		boxz = atof(box_z.c_str());
	    }
	    //
	}	
    }
    OpenFile.close();
    int t ;
    t = ii - 1 ;
    cout<<"####"<<t<<"####"<< Xu[t] << "," << typeXu[t] << "," << x[t] << "," << y[t] << "," << z[t] <<endl ;
    cout<<boxx<<":"<<endl;
    cout<<boxy<<":"<<endl;
    cout<<boxz<<":"<<endl;
    double const PI = acos(double(-1));
    cout<<PI<<endl;
    //read_size
    ifstream OpenFileB("DHOP35du.txt", ios::in);
    if (OpenFileB.fail())
    {
        cout<<"Can not open target file"<<endl;
        return 0;
    }
    string csstr;
    int crystal_size = 0;
    std::vector<int> crystal_size_index;
    int coutttt = 0;
    int coutttt_B = 0;
    int coutttt_A = 0;
    double centre_x = 0;
    double centre_y = 0;
    double centre_z = 0;
    while (OpenFileB >> csstr)
    {
	coutttt ++;
	if (coutttt == 1) {
	   coutttt_A = atoi(csstr.c_str());
	   cout<<coutttt_A<<endl; 
	}
	if ((coutttt >2) && (coutttt<=(coutttt_A+2)) ) {
	    crystal_size ++;
	    crystal_size_index.push_back(atoi(csstr.c_str())) ;
//	    cout<<coutttt<<" | "<< coutttt_B<<" | "<< csstr <<" "<<crystal_size_index[coutttt_B]<<" x: "<<x[crystal_size_index[coutttt_B]]<<" y: "<<y[crystal_size_index[coutttt_B]]<<" z: "<<z[crystal_size_index[coutttt_B]]<<endl;
	    centre_x = centre_x + x[crystal_size_index[coutttt_B]] ;
	    centre_y = centre_y + y[crystal_size_index[coutttt_B]] ;
	    centre_z = centre_z + z[crystal_size_index[coutttt_B]] ;  
	    coutttt_B ++ ;
	}
    } 
    centre_x = centre_x / coutttt_B ;
    centre_y = centre_y / coutttt_B ; 
    centre_z = centre_z / coutttt_B ; 
    cout<<coutttt_B<<" "<<centre_x<<" "<<centre_y<<" "<<centre_z<<endl;
    OpenFileB.close();
    ////////////////////////////////////////////////////////////////////////////////////////
    double distance_BF = 0 ;
    double distance_maxF = 0;
    double distance_minF = 10000;
    int d_min_index = 100000;
    int d_max_index = 100000;
    double crystal_size_index_x[coutttt_B] ;
    double crystal_size_index_y[coutttt_B] ;
    double crystal_size_index_z[coutttt_B] ;
    int crystal_size_indexF[coutttt_B];
    for ( int i_max = 0 ; i_max < coutttt_B ; i_max ++ ){
    	crystal_size_indexF[i_max] = crystal_size_index[i_max];
    	crystal_size_index_x[i_max] = x[crystal_size_index[i_max]] ;
    	crystal_size_index_y[i_max] = y[crystal_size_index[i_max]] ;
    	crystal_size_index_z[i_max] = z[crystal_size_index[i_max]] ;
	distance_BF = ( ( centre_x - x[crystal_size_index[i_max]] ) * ( centre_x - x[crystal_size_index[i_max]] ) ) + ( ( centre_y - y[crystal_size_index[i_max]] ) * ( centre_y - y[crystal_size_index[i_max]] ) ) + ( ( centre_z - z[crystal_size_index[i_max]] ) * ( centre_z - z[crystal_size_index[i_max]] ) ) ;
	if ( distance_minF > distance_BF ) {
	   distance_minF = distance_BF ;
	   d_min_index = crystal_size_index[i_max];
	}
	if ( distance_maxF < distance_BF ) {
	   distance_maxF = distance_BF ;
	   d_max_index = crystal_size_index[i_max];
	}
    } 
    cout<<"min: "<<distance_minF<<" "<<sqrt(distance_minF)<<" ["<<d_min_index<<"] "<<"max: "<<distance_maxF<<" "<<sqrt(distance_maxF)<<" ["<<d_max_index<<"] "<<endl;
   // cout<<"===============minx,y,z"<<crystal_size_index_x[d_min_index]<<"="<<crystal_size_index_y[d_min_index]<<"="<<crystal_size_index_z[d_min_index]<<endl;
    cout<<"======================="<<endl;
    CENTRE (coutttt_B , &centre_x , &centre_y , &centre_z , crystal_size_index_x,crystal_size_index_y,crystal_size_index_z );
    cout<<"&centre_x "<<centre_x<<"&centre_y "<<centre_y<<"&centre_z "<<centre_z<<endl;
    //cout<<"===============minx,y,z"<<crystal_size_index_x[d_min_index]<<"="<<crystal_size_index_y[d_min_index]<<"="<<crystal_size_index_z[d_min_index]<<endl;
    cout<<"======================="<<endl;
    //cout<<"===============minx,y,z"<<crystal_size_index_x[d_min_index]<<"="<<crystal_size_index_y[d_min_index]<<"="<<crystal_size_index_z[d_min_index]<<endl;
    //
    MIN_MAX_CENTRE(&distance_BF,&distance_maxF,&distance_minF,&d_min_index,&d_max_index,coutttt_B,&centre_x,&centre_y,&centre_z,crystal_size_indexF,crystal_size_index_x,crystal_size_index_y,crystal_size_index_z);
    cout<<"===============minx,y,z"<<crystal_size_index_x[d_min_index]<<"="<<crystal_size_index_y[d_min_index]<<"="<<crystal_size_index_z[d_min_index]<<endl;
    //
    if ( sqrt(distance_minF) < 0.35 ) {
	//
	cout<<abs(crystal_size_index_x[d_max_index] - centre_x)<<endl;
	cout<<abs(crystal_size_index_y[d_max_index] - centre_y)<<endl; 
	cout<<abs(crystal_size_index_z[d_max_index] - centre_z)<<endl; 
	if ( ( abs(crystal_size_index_x[d_max_index] - centre_x ) < (boxx/2) ) && ( abs(crystal_size_index_y[d_max_index] - centre_y) < (boxy/2) ) && ( abs(crystal_size_index_z[d_max_index] - centre_z ) < (boxz/2) ) ) {
	}
	else{
	//for()_>new_centre_X
	     for ( int re_i = 0 ; re_i < coutttt_B ; re_i ++ ) {
		  double linshi_x = crystal_size_index_x[re_i] - centre_x;
		  double linshi_y = crystal_size_index_y[re_i] - centre_y;
		  double linshi_z = crystal_size_index_z[re_i] - centre_z;
		  if ( linshi_x > (boxx/2) ) {
			  crystal_size_index_x[re_i] = crystal_size_index_x[re_i] - boxx ;
		  }
		  if ( linshi_x < (-1 * boxx/2) ) {
			  crystal_size_index_x[re_i] = crystal_size_index_x[re_i] + boxx ;
		  }
		  ////////////////////////////////////////_________________________________________________x
		  if ( linshi_y > (boxy/2) ) {
			  crystal_size_index_y[re_i] = crystal_size_index_y[re_i] - boxy ;
		  }
		  if ( linshi_y < (-1 * boxy/2) ) {
			  crystal_size_index_y[re_i] = crystal_size_index_y[re_i] + boxy ;
		  }
		  ///////////////////////////////////////__________________________________________________y
		  if ( linshi_z > (boxz/2) ){
			  crystal_size_index_z[re_i] = crystal_size_index_z[re_i] - boxz ;
		  }
		  if ( linshi_z < (-1 * boxz/2) ){
			  crystal_size_index_z[re_i] = crystal_size_index_z[re_i] + boxz ;
		  }
		  //////////////////////////////////////////
	     }
	     //
	     cout<<"======================="<<endl;
	     CENTRE (coutttt_B , &centre_x , &centre_y , &centre_z , crystal_size_index_x,crystal_size_index_y,crystal_size_index_z );
             cout<<"&centre_x "<<centre_x<<"&centre_y "<<centre_y<<"&centre_z "<<centre_z<<endl;
	     cout<<"======================="<<endl;
	     MIN_MAX_CENTRE(&distance_BF,&distance_maxF,&distance_minF,&d_min_index,&d_max_index,coutttt_B,&centre_x,&centre_y,&centre_z,crystal_size_indexF,crystal_size_index_x,crystal_size_index_y,crystal_size_index_z);
	     cout<<": "<<crystal_size_index_x[d_max_index]<<" ; "<<crystal_size_index_y[d_max_index]<<" ; "<<crystal_size_index_z[d_max_index]<<endl;
	     //
	}
    }
    else{
    	centre_x = crystal_size_index_x[d_min_index] ;
	centre_y = crystal_size_index_y[d_min_index] ; 
     	centre_z = crystal_size_index_z[d_min_index] ;
	cout<<" else:new x,yz :"<<centre_x<<" , "<<centre_y<<" , "<<centre_z<<endl;
        //for()_>new_ce
	for( int re_i = 0 ; re_i < coutttt_B ; re_i ++ ) {
	     double linshi_x = crystal_size_index_x[re_i] - centre_x;
             double linshi_y = crystal_size_index_y[re_i] - centre_y;
             double linshi_z = crystal_size_index_z[re_i] - centre_z;
	     if ( linshi_x > (boxx/2) ) {
	         crystal_size_index_x[re_i] = crystal_size_index_x[re_i] - boxx ;
	     }
	     if ( linshi_x < (-1 * boxx/2) ) {
	         crystal_size_index_x[re_i] = crystal_size_index_x[re_i] + boxx ;
             }
            ////////////////////////////////////////_______________________________________________
	     if ( linshi_y > (boxy/2) ) {
                 crystal_size_index_y[re_i] = crystal_size_index_y[re_i] - boxy ;
             }
             if ( linshi_y < (-1 * boxy/2) ) {
                 crystal_size_index_y[re_i] = crystal_size_index_y[re_i] + boxy ;
             }
             ///////////////////////////////////////__________________________________________________y
             if ( linshi_z > (boxz/2) ){
	         crystal_size_index_z[re_i] = crystal_size_index_z[re_i] - boxz ;
	     }
	     if ( linshi_z < (-1 * boxz/2) ){
	         crystal_size_index_z[re_i] = crystal_size_index_z[re_i] + boxz ;
	     }
	     //////////////////////////////////////////
	}
	cout<<"else======================="<<endl;
	CENTRE (coutttt_B , &centre_x , &centre_y , &centre_z , crystal_size_index_x,crystal_size_index_y,crystal_size_index_z );
	cout<<"&centre_x "<<centre_x<<"&centre_y "<<centre_y<<"&centre_z "<<centre_z<<endl;
        cout<<"else======================="<<endl;
        MIN_MAX_CENTRE(&distance_BF,&distance_maxF,&distance_minF,&d_min_index,&d_max_index,coutttt_B,&centre_x,&centre_y,&centre_z,crystal_size_indexF,crystal_size_index_x,crystal_size_index_y,crystal_size_index_z);
        cout<<": "<<crystal_size_index_x[d_max_index]<<" ; "<<crystal_size_index_y[d_max_index]<<" ; "<<crystal_size_index_z[d_max_index]<<"==="<<endl;
	//if()
	if ( ( abs(crystal_size_index_x[d_max_index] - centre_x ) < (boxx/2) ) && ( abs(crystal_size_index_y[d_max_index] - centre_y) < (boxy/2) ) && ( abs(crystal_size_index_z[d_max_index] - centre_z ) < (boxz/2) ) ) {
	}
	else{
	     for ( int re_i = 0 ; re_i < coutttt_B ; re_i ++ ) {
	  	  double linshi_x = crystal_size_index_x[re_i] - centre_x;
		  double linshi_y = crystal_size_index_y[re_i] - centre_y;
		  double linshi_z = crystal_size_index_z[re_i] - centre_z;
		  if ( linshi_x > (boxx/2) ) {
	    	       crystal_size_index_x[re_i] = crystal_size_index_x[re_i] - boxx ;
		  }
		  if ( linshi_x < (-1 * boxx/2) ) {
		       crystal_size_index_x[re_i] = crystal_size_index_x[re_i] + boxx ;
		  }
		   ////////////////////////////////////////_________________________________________________x
		  if ( linshi_y > (boxy/2) ) {
		       crystal_size_index_y[re_i] = crystal_size_index_y[re_i] - boxy ;
		  }
		  if ( linshi_y < (-1 * boxy/2) ) {
		       crystal_size_index_y[re_i] = crystal_size_index_y[re_i] + boxy ;
		  }
		  ///////////////////////////////////////__________________________________________________y
		  if ( linshi_z > (boxz/2) ){
		        crystal_size_index_z[re_i] = crystal_size_index_z[re_i] - boxz ;
		  }
		  if ( linshi_z < (-1 * boxz/2) ){
		        crystal_size_index_z[re_i] = crystal_size_index_z[re_i] + boxz ;
		  }
		  //////////////////////////////////////////
	     }
	     //
	     cout<<"new======================="<<endl;
	     CENTRE (coutttt_B , &centre_x , &centre_y , &centre_z , crystal_size_index_x,crystal_size_index_y,crystal_size_index_z );
	     cout<<"&centre_x "<<centre_x<<"&centre_y "<<centre_y<<"&centre_z "<<centre_z<<endl;
	     cout<<"nwe======================="<<endl;
	     MIN_MAX_CENTRE(&distance_BF,&distance_maxF,&distance_minF,&d_min_index,&d_max_index,coutttt_B,&centre_x,&centre_y,&centre_z,crystal_size_indexF,crystal_size_index_x,crystal_size_index_y,crystal_size_index_z);
	     cout<<": "<<crystal_size_index_x[d_max_index]<<" ; "<<crystal_size_index_y[d_max_index]<<" ; "<<crystal_size_index_z[d_max_index]<<endl;
	     //
	}	     
    }
    double hf = 0 -91.093 ;
    cout << abs(hf) <<endl;
    //
    Eigen::Vector3d v_3d_cs;
    Eigen::Matrix3d matrix_33_cs = Eigen::Matrix3d::Zero();
    //
   // v_3d_cs << centre_x - 1,centre_y - 1,centre_z - 1 ;
   // matrix_33_cs = v_3d_cs * v_3d_cs.transpose() ;
   // cout<<"v_3d_cs:\n"<<v_3d_cs<<endl;
   // cout<<"matrix_33_cs:\n"<<matrix_33_cs<<endl;
  //  Eigen::Matrix3d matrix_33_csB = Eigen::Matrix3d::Zero();
  //  matrix_33_csB << 1,1,1,1,1,1,2,2,2;
  //  matrix_33_cs = matrix_33_cs + matrix_33_csB ;
   // cout<<"matrix_33_cs:\n"<<matrix_33_cs<<endl;
    for ( int re_i = 0 ; re_i < coutttt_B ; re_i ++ ) {
    	v_3d_cs << crystal_size_index_x[re_i] - centre_x , crystal_size_index_y[re_i] - centre_y, crystal_size_index_z[re_i] - centre_z ;
//	cout<<"v_3d_cs:\n"<<v_3d_cs<<endl;
//	cout<<"v_3d_cs*T:\n"<<v_3d_cs* v_3d_cs.transpose() <<endl;
	matrix_33_cs = matrix_33_cs + ( v_3d_cs * v_3d_cs.transpose() );
    }
    cout<<"matrix_33_cs:\n"<<matrix_33_cs<<endl;
    matrix_33_cs = matrix_33_cs / coutttt_B ;
    cout<<"matrix_33_cs:\n"<<matrix_33_cs<<endl;
    //for (  ) {
    //}
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>eigen_solver(matrix_33_cs.transpose()*matrix_33_cs);
    cout << "Eigen value = \n" << eigen_solver.eigenvalues() << endl;
    cout << "Eigen vectors = \n" << eigen_solver.eigenvectors() << endl;
    Eigen::Vector3d v_3d_fin;
    v_3d_fin = eigen_solver.eigenvalues() ;
    cout << " v_3d_fin \n"<<v_3d_fin(0,0)<<" "<<v_3d_fin(1,0)<<" "<<v_3d_fin(2,0)<<endl;
    cout<<"================================"<<endl;
    EigenSolver<Matrix3d> ePPP(matrix_33_cs);
    Matrix3d DPPP = ePPP.pseudoEigenvalueMatrix();
    cout<<"DPPP:\n"<<DPPP<<endl;
    cout<<"================================"<<endl;
    ofstream OUT_print_QD;
    OUT_print_QD.open("asphericity.txt",ios::trunc);
    double F_qiu ;
    double F_a = ((pow(v_3d_fin(0,0),2))+(pow(v_3d_fin(1,0),2))+(pow(v_3d_fin(2,0),2))) ;
    double F_b = ((v_3d_fin(0,0)+v_3d_fin(1,0)+v_3d_fin(2,0))*(v_3d_fin(0,0)+v_3d_fin(1,0)+v_3d_fin(2,0))) ;
    cout<<((pow(v_3d_fin(0,0),2))+(pow(v_3d_fin(1,0),2))+(pow(v_3d_fin(2,0),2)))<<endl;
    cout<<((v_3d_fin(0,0)+v_3d_fin(1,0)+v_3d_fin(2,0))*(v_3d_fin(0,0)+v_3d_fin(1,0)+v_3d_fin(2,0)))<<endl;
    cout<<(3.0/2.0)*(F_a/F_b)<<endl;
    F_qiu = ( (3.0/2.0) * ( ((pow(v_3d_fin(0,0),2))+(pow(v_3d_fin(1,0),2))+(pow(v_3d_fin(2,0),2))) / ((v_3d_fin(0,0)+v_3d_fin(1,0)+v_3d_fin(2,0))*(v_3d_fin(0,0)+v_3d_fin(1,0)+v_3d_fin(2,0))) ) ) - 0.5;
    cout<<"F_qiu : "<<F_qiu <<endl;
    cout<<1 +(2*2* boxx/2 ) <<endl;
    OUT_print_QD<<F_qiu;
    OUT_print_QD.close();
   // double PO = pow(2,3);
  //  cout<<PO<<endl;

    return 0;
}
