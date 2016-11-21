/*
 *  Name: Xiao Chen
 *  Stony Brook ID:111156560
 *  Final Project First Week Progress
 *  This program included the function of reading three files and all the equations that needed to be used.
 *  However, the calculation has not been done yet.
 *  The test of reading file and putting the data into the three turbine arrays has finished and would give a output of the x, y, a and row for all the turbines.
 *  In this version, the power can be calculated and output.
 *  However, there are still some bugs there.
 *  The results output might not correct.
 *  The bugs might be correct in the next version.
 */
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

#define pi 3.14159
#define k 0.4
#define f 1e-4
#define A 4
#define density 1.225
using namespace std;

int num_of_turbine;

class Turbine{//turbine properties read from file
protected:
    int horns_rev;
    int wind_tunnel;
    double bl_height;
    double c_alfa;//eq.13
    int num_of_turbine;
    double geo_wind_speed;
    double ground_roughness_height;
    double diameter;
    double tur_hub_height;
    double sx;
    double sy;
    double fraction;//eq.7
    double location_ibl_start;
    double beta;//eq.10
    double x;
    double y;
    double a;
    int row;
    int turbine_num;
public:
    void set_num(int num){
        turbine_num = num;
    }
    void set_x(double x1){
        x = x1;
    }
    void set_y(double y1){
        y = y1;
    }
    void set_a(double a1){
        a = a1;
    }
    void set_row(double row1){
        row = row1;
    }
    double get_row(){
        return row;
    }
    double get_x(){
        return x;
    }
    double get_y(){
        return y;
    }
    void show_x(){
        cout << x << endl;
    }
    void set_parameter(double parameter[]){
        horns_rev = parameter[0];
        wind_tunnel = parameter [1];
        bl_height = parameter [2];
        c_alfa = parameter [3];
        num_of_turbine = parameter [4];
        geo_wind_speed = parameter [5];
        ground_roughness_height =parameter [6];
        diameter = parameter [7];
        tur_hub_height =parameter[8];
        sx = parameter[9];
        sy = parameter[10];
        fraction = parameter [11];
        location_ibl_start= parameter[12];
        beta = parameter[13];
            //para_arr[i] = parameter[i];
        //cout << tur_hub_height << endl;
    }
    void show(){
        cout << x << "\t" << y << "\t" << a << "\t" <<row << endl;
    }
};

class Turbine_dev : public Turbine{//turbine properties needed to be calculated
protected:
    double ct;
    double cp;
    double cft;
    double z0hi;
    double u_0;
    double u_hi;
    double hg;
    double xfd;
    double uh;
    double uh_fd;
    double ua_fd;
    double incoming_wind_speed;
    double u_a ;// u*a
    double ua;//Ua
    double u_nw;
    double x_for_cal;
    double rw;
    double alpha;
    double downwind_dis;
    double velocity_in_wake;
    double power;
    double du_sum ;
    double velocity;
    double du;
    double old_velocity;

public:
    void compute_ct(){// compute ct 1
        ct = 4*a*(1-a);
        //cout << "ct " << ct << endl;
    }
    void compute_cp(){// compute cp 1
        cp = 4*a*pow((1-a),2);
       // cout << "cp " << cp << endl;
    }
    void compute_cft(){//compute cft and use ct 2
        cft = ct*pi*pow(diameter,2)/(4*sx*sy);
        //cout << "cft " << cft << endl;
    }
    void compute_z0hi(){//compute z0hi and use cft 3
        z0hi = tur_hub_height*exp(-k/sqrt(cft/2+pow(k/log(tur_hub_height/ground_roughness_height),2)));
        //cout << "z0hi " << z0hi << endl;
    }
    void compute_u_0(){//compute u*0 1
        if(turbine_num == 3){
            u_0 = k*geo_wind_speed/(log(geo_wind_speed/(f*exp(A)*ground_roughness_height)));
        }else{
            u_0 = k*geo_wind_speed/(log(bl_height/ground_roughness_height));
        }
       // cout << "u_0 " << u_0 << endl;

    }
    void compute_u_hi(){//compute u*hi and use z0hi 4
        if(turbine_num == 3){
            u_hi = k*geo_wind_speed/(log(geo_wind_speed/(f*exp(A)*z0hi)));
        }else{
            u_hi = k*geo_wind_speed/(log(bl_height/z0hi));
        }
    }
    void compute_hg(){// compute hg and use u*0, u*hi 5
        hg = exp((u_hi*log(z0hi)-u_0*log(ground_roughness_height))/(u_hi-u_0));
        //cout << "hg " << hg << endl;
    }
    void compute_xfd(){//!!compute xfd and use hg 6 (may be no correct)
        xfd = pow(hg-tur_hub_height,5/4)*pow(z0hi,-1/4);

        //cout << "xfd " << xfd <<endl;
    }
    void compute_uh(){// compute uh and use u*0 2
        uh =u_0*log(tur_hub_height/ground_roughness_height)/k;
       // cout << "uh "<<uh << endl;
       // cout <<"tub heigh "<< tur_hub_height<<endl;
       // cout << "ground roughness "<< ground_roughness_height<<endl;
       // cout << "k " << k << endl;
    }
    void compute_uh_fd(){//compute uh_fd and use u*hi 5
        uh_fd = u_hi/k*log(tur_hub_height/z0hi);
    }
    void compute_ua_fd(){// compute ua_fd and use uh, uhfd 6
        double gamma = fraction;
        ua_fd = gamma*uh_fd + (1-gamma)*uh;
    }
    void initialize_old_velocity(){//initialize incoming wind speed use uh 3
        old_velocity = uh;
    }

    void initialize_all_value(){
        compute_ct();//1
        compute_cp();//1
        compute_u_0();//1
        compute_cft();//2
        compute_uh();// 2
        compute_z0hi();//3
        initialize_old_velocity();//3
        compute_u_hi();//4
        compute_hg();//5
        compute_uh_fd();//5
        compute_xfd();//6
        compute_ua_fd();//6
        initial_first_row();
    }
    void initial_first_row(){
        if(row == 1){
            ua = uh;
            velocity = uh;
            incoming_wind_speed = uh;
            u_nw = (1-2*a)*incoming_wind_speed;
            compute_u_a();
            alpha = c_alfa*u_a/u_nw;

        }
    }
    //=========================================================================
    void compute_x(Turbine_dev T){//compute x for calculate 1
        x_for_cal = abs(x-T.x);
        //cout << row << "\t "<< x_for_cal << endl;
    }

    void compute_downwind_dis(Turbine_dev T){//compute downwind distance 1
        downwind_dis = abs(x - T.get_x());
        //cout << row << "\t"<<downwind_dis << endl;
    }

    void compute_alpha(){//compute alpha(use u_a, u_nw)//8
        alpha = c_alfa*u_a/u_nw;
        //cout << "u_a" << u_a <<" ";
    }
    void compute_unw(){//compute unw //3
        if(row == 1){
            u_nw = (1-2*a)*uh;
        }else{
            u_nw = (1-2*a)*incoming_wind_speed;
        }
    }
    void compute_u_a(){//compute u_a (use xfd, x for cal)//7
        if( x_for_cal >= 0 && x_for_cal < beta*xfd){
            u_a = u_0;
        }else if(x_for_cal>= beta*xfd && x_for_cal < xfd){
            u_a = (1-(x_for_cal-beta*xfd)/((1-beta)*xfd)*u_0 + (x_for_cal-beta*xfd)/((1-beta)*xfd)*u_hi);
        }else{
            u_a = u_hi;
        }
    }
    void compute_power(){//compute power
        if(row == 1){
            power = cp*density*pi*pow(diameter/2,2)*pow(incoming_wind_speed,3)/2;
           // cout << "cp" <<cp << endl;
        }else{
            power = cp*density*pi*pow(diameter/2,2)*pow(velocity,3)/2;
        }     //cout << incoming_wind_speed <<" "<<velocity_in_wake << " " << alpha<<endl;
        //cout << rw <<" " << x_for_cal<< " " << xfd << " "<< beta*xfd <<endl;
    }

    void compute_du(Turbine_dev T){ //compute du // 10
        if(abs(y-T.y)>rw){
            du = 0;
        }else{
            du = incoming_wind_speed - old_velocity;
        }
        velocity = old_velocity + du;
        old_velocity = velocity;
    }
    void compute_ua(){//compute Ua, use ua_fd 7
        if(x_for_cal>=0 && x_for_cal < xfd){
            ua = (1-x_for_cal/xfd)*uh+x_for_cal/xfd*ua_fd;
        }else{
            ua = ua_fd;
        }
    }
    void compute_rw(){//because of alpha (use alpha)// 9
        if(downwind_dis < diameter){
            rw = diameter/2;
        }else{
            rw = diameter/2+alpha*(downwind_dis-diameter);
        }
       // cout <<row << "\t" <<rw << " " << endl;
      // cout << row << " " << alpha << endl;
    }
    void compute_velocity_in_wake(){//(use ua, rw) for this turbine // 10
        if(downwind_dis < diameter){
            velocity_in_wake = u_nw;
        }else{
            velocity_in_wake = (u_nw + ua*sqrt(pow(u_nw/ua,2)+4*(1-pow((diameter/2)/rw,2))*(1-u_nw/ua)))/2;
        }
       //cout << row <<"\t"<<velocity_in_wake << "\t" << rw<< endl;
    }
    void compute_incoming_velocity(Turbine_dev T){
        incoming_wind_speed = T.velocity_in_wake;
    }
    void show_power(){
        cout << "Power in "<<row << " is " <<power << endl;
        //cout << row <<" x " << x_for_cal << " xfd " << xfd << " ua " << ua << endl;
    }
};

void readfile(Turbine_dev turbine_array[], int num); // read file
void set_value_xyarow(double x[], double y[], double a[], double row[],Turbine_dev T[],int num2,int filenum);// set value of x, y, a and row of turbine
void read_turbine_num(int num,int &turbine_num);// read the number of turbine from file
void compute_power(Turbine_dev *T);
void sort_turb(Turbine_dev *T);

int main(){
    int type;
    char yesorno;
    do{
        cout << "Please input Type (1-Align, 2-Stagger, 3-Horn Rev):";
        cin >> type;
        if(type == 1){
            cout << "Power of WindTunnel-Align: "<< endl;
            read_turbine_num(1,num_of_turbine);// get the number of turbines from file 1
            Turbine_dev *T1;
            T1 = new Turbine_dev [num_of_turbine];
            readfile(T1,1);

        //------------------------------------------------------------------------
        compute_power(T1);
        delete []T1;
        }else if(type == 2){
            cout << "Power of WindTunnel-Stagger: "<< endl;
        //====================================================================
            read_turbine_num(2,num_of_turbine);// get the number of turbines from file 1

            Turbine_dev *T2;
            T2 = new Turbine_dev [num_of_turbine];

            readfile(T2,2);
            compute_power(T2);
            delete []T2;
        }else if(type ==3){
            cout << "Power of Horn Rev: "<< endl;
        //--------------------------------------------------------------
            read_turbine_num(3,num_of_turbine);// get the number of turbines form file 3
            Turbine_dev *T3;
            T3 = new Turbine_dev [num_of_turbine];
            readfile(T3,3);
            compute_power(T3);
            delete []T3;
        }
    //--------------------------------------------------------------
        cout << "\nTry again?(y/n)";
        cin >> yesorno;
    }while(yesorno=='y');
    return 0;
}
void read_turbine_num(int num, int &turbine_num){//read the turbine number from file
        string str;
        string temp;
        double tub_num;
        switch(num){
            case 1:
                str = "WindTunnel_Align.txt";
                break;
            case 2:
                str = "WindTunnel_Stagger.txt";
                break;
            case 3:
                str = "HornsRev.txt";
                break;
        }
        ifstream myReadFile;

        myReadFile.open(str.c_str());
        if(myReadFile.is_open()){
            while(!myReadFile.eof()){
                getline(myReadFile,temp);
                if(temp.compare("Number of Turbines")== 0){
                    myReadFile >> turbine_num;
                }
            }
        }
}

void readfile(Turbine_dev turbine_array[], int num){//read parameter from file
        string str;
        string temp;
        int lines = 0;
        switch(num){
            case 1:
                str = "WindTunnel_Align.txt";
                break;
            case 2:
                str = "WindTunnel_Stagger.txt";
                break;
            case 3:
                str = "HornsRev.txt";
                break;
        }
        ifstream myReadFile;

        myReadFile.open(str.c_str());

        if(myReadFile.is_open()){
            while(!myReadFile.eof()){
                getline(myReadFile,temp);
                lines++;
            }
        }else{
            cout << "Cannot Open File." << endl;
        }
        myReadFile.close();

        myReadFile.open(str.c_str());

        if(myReadFile.is_open()){
            double *parameter;
            double *x_p;
            double *y_p;
            double *a_p;
            double *row_p;

            parameter = new double[lines * 10];
            x_p = new double [lines];
            y_p = new double [lines];
            a_p = new double [lines];
            row_p = new double [lines];
            int current_line = 1;

            bool startTubineInformation = false;
            int num1 = 0;
            int num2 = 0;
            while(current_line != lines) {
                      if(!startTubineInformation)
                      {
                        getline(myReadFile, temp);
                        current_line++;
                        //cout << temp << endl;
                        if(temp.compare("Tubine information") == 0)
                        {
                            startTubineInformation = true;
                            getline(myReadFile, temp);   // x	y	a	row#
                            current_line++;
                           // cout << temp << endl;
                        }
                      }

                        // cout << startTubineInformation << endl;
                      if(startTubineInformation)
                      {
                            myReadFile >> x_p[num2]>> y_p[num2] >> a_p[num2] >> row_p[num2]; // x  y  z row#
                           // cout << x_p[num2]<< "\t"<<y_p[num2] <<"\t"<< a_p[num2] <<"\t"<< row_p[num2] << endl;
                            num2 ++;
                      }
                      else
                      {
                          myReadFile >> parameter[num1];
                         // cout << parameter[num1] << endl;
                          num1 ++;
                      }

                      for(int i = 0; i < num_of_turbine; i++){
                        turbine_array[i].set_parameter(parameter);
                      }

                      getline(myReadFile,temp);
                      current_line++;
                }
                set_value_xyarow(x_p,y_p,a_p,row_p,turbine_array,num2,num);
                delete []x_p;
                delete []y_p;
                delete []a_p;
                delete []row_p;

        }else{
            cout << "Cannot Open File." << endl;
        }

        myReadFile.close();
    }

void set_value_xyarow(double x[], double y[], double a[], double row[],Turbine_dev T[],int num2,int filenum){// set value of x, y, a ,row
    for (int i=0; i < num2; i ++){
        T[i].set_x(x[i]);
        T[i].set_y(y[i]);
        T[i].set_a(a[i]);
        T[i].set_row(row[i]);
        T[i].set_num(filenum);
    }
}

void compute_power(Turbine_dev *T){

        int NumofTurbineinRow []= {0,0,0,0,0,0,0,0,0,0};
        int row;
        for(int i=0; i < num_of_turbine; i++){
            row = (int)T[i].get_row();
            NumofTurbineinRow[row-1]++;
        }
        sort_turb(T);
        int orign_row = T[0].get_row();

        for(int i =0; i <num_of_turbine; i ++){
            T[i].initialize_all_value();
        }// initialize all value

        for(int i =0; i < num_of_turbine; i ++){
            for(int j = 0; j < num_of_turbine;j++){
                if(T[j].get_row() == 1){
                    T[i].compute_x(T[j]);
                }
            }
        }// initialize all x_for_cal
        for(int i = 0 ; i < num_of_turbine; i++){
            T[i].initial_first_row();
        }// initialize first row


        for(int i = 0 ; i < num_of_turbine; i++){
            for(int j = 0; j< num_of_turbine ; j++){
                if(T[i].get_row()<T[j].get_row()){
                    T[j].compute_downwind_dis(T[i]);
                    T[i].compute_downwind_dis(T[j]);
                    T[i].compute_u_a();
                    T[i].compute_unw();
                    T[i].compute_alpha();
                    T[i].compute_rw();
                    T[i].compute_ua();
                    T[i].compute_velocity_in_wake();
                    T[j].compute_incoming_velocity(T[i]);
                    T[j].compute_u_a();
                    T[j].compute_unw();
                    T[j].compute_alpha();
                    T[j].compute_rw();
                    T[j].compute_ua();
                    T[j].compute_du(T[i]);
                }
            }
        }

        for(int i = 0; i < num_of_turbine; i ++){
                T[i].compute_power();
        }
        for(int i =0; i< num_of_turbine; i++){
            T[i].show_power();
        }
}
void sort_turb(Turbine_dev *T){

    Turbine_dev temp_turb;
    for(int i=0; i < num_of_turbine; i ++){
        for(int j =i; j < num_of_turbine; j++){
            if(T[i].get_row()> T[j].get_row()){
                temp_turb = T[i];
                T[i] = T[j];
                T[j] = temp_turb;
            }
        }
    }

}
