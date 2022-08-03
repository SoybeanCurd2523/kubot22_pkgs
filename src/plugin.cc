#include <gazebo/common/common.hh>
#include <gazebo/common/Plugin.hh>
#include <ros/ros.h>
#include <boost/bind.hpp>
#include <gazebo/gazebo.hh>
#include <gazebo/physics/physics.hh>
#include <gazebo/sensors/sensors.hh>
#include <stdio.h>
#include <iostream>
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/String.h>
#include <functional>
#include <ignition/math/Vector3.hh>

#include "Eigen/Dense"
#include <rbdl/rbdl.h>
#include <rbdl/addons/urdfreader/urdfreader.h>

#include <fstream>

#include "gnuplot.h"

#define PI 3.141592
#define D2R 3.141592/180
#define R2D 180/PI


using Eigen::MatrixXd;
using Eigen::VectorXd;

using Eigen::Vector3d;
using Eigen::RowVectorXd;
using Eigen::RowVector3d;

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace std;

Model* kubot22_Model = new Model;
VectorXd calG(12);
VectorXd calC(12);


RigidBodyDynamics::Math::VectorNd RobotState(12);
RigidBodyDynamics::Math::VectorNd RobotStatedot(12);
RigidBodyDynamics::Math::VectorNd RobotStatedot2(12);


VectorXd Tau = VectorXd::Zero(12);
VectorXd Nowq = VectorXd::Zero(12);
VectorXd NowDq = VectorXd::Zero(12);
VectorXd NowD2q = VectorXd::Zero(12);


float LJ[6] = {0};
float RJ[6] = {0};
 
float N_LJ[6] = {0};
float N_RJ[6] = {0};
 
int CNT = 0;
double T = 3;

VectorXd q_init(6), q_cal(6), q0(6), q(6), q0_init(6), q0_final(6), q_prev(6);
VectorXd q_init_R(6), q_cal_R(6), q0_R(6), q_R(6), q0_init_R(6), q0_final_R(6), q_prev_R(6);


double cnt = 0;



///////////// 선언 ////////////

double All_time_trajectory;
double mydt;
const int N = 1000;
int n;

double z_c;
double g ;

MatrixXd A(3,3);
VectorXd B(3);
RowVectorXd C(3);

double Gi;

void getGi(){
    FILE *fp;
    fp = fopen("/home/ubuntu/catkin_ws/src/kubot22_pkgs/src/Gi.txt", "r"); // 파일 열기

    if(fp == NULL){
        ROS_INFO("Gi read failed");
        return;
    }

    fscanf(fp, "%lf", &Gi);
    fclose(fp);

    ROS_INFO("read Gi value : %lf", Gi);    
}

RowVectorXd Gx(3);

double Gx_1, Gx_2, Gx_3;

void getGx(){
    FILE *fp;
    fp = fopen("/home/ubuntu/catkin_ws/src/kubot22_pkgs/src/Gx.txt", "r"); // 파일 열기

    if(fp == NULL){
        ROS_INFO("Gx read failed");
        return;
    }

    fscanf(fp, "%lf %lf %lf", &Gx_1, &Gx_2, &Gx_3);

    fclose(fp);

    ROS_INFO("read Gx[0] value : %lf", Gx_1);    
}


double Gp[N];

void getGp(){
    FILE *fp;
    fp = fopen("/home/ubuntu/catkin_ws/src/kubot22_pkgs/src/Gp.txt", "r"); // 파일 열기

    if(fp == NULL){
        ROS_INFO("Gp read failed");
        return;
    }

    for(int i=0 ; i<N ; i++){
        fscanf(fp, "%lf", &Gp[i]);
    }
    fclose(fp);

    ROS_INFO("read Gp[0] value : %lf", Gp[0]);    
}


VectorXd X(3);
VectorXd X_new(3);
RowVectorXd zmp_x_ref(n); // 1*15001
double u_x;

VectorXd Y(3);
VectorXd Y_new(3);
RowVectorXd zmp_y_ref(n);
double u_y;

RowVectorXd zmp_x(n);
double zmp_x_old;

RowVectorXd zmp_y(n);
double zmp_y_old = 0;

RowVectorXd com_x(n);

RowVectorXd com_y(n);

double err_x;
double err_y;

double sum_e_x;
double sum_e_y;

double u_sum_p_x;
double u_sum_p_y;

double mytime;
int i=0;

ofstream fout;


void setCom(int k){
    X = X_new;
    Y = Y_new;

    err_x = zmp_x_old - zmp_x_ref(k);
    err_y = zmp_y_old - zmp_y_ref(k);

    sum_e_x = sum_e_x + err_x;
    sum_e_y = sum_e_y + err_y;

    for(int j=0 ; j< N ; j++){
        u_sum_p_x = u_sum_p_x + Gp[j] * zmp_x_ref(k+j);
        u_sum_p_y = u_sum_p_x + Gp[j] * zmp_y_ref(k+j -0.001);
    }

    u_x = -Gi * sum_e_x - Gx * X - u_sum_p_x;
    u_y = -Gi * sum_e_y - Gx * Y - u_sum_p_y;

    X_new = A*X + B*u_x;
    Y_new = A*Y + B*u_y;

    zmp_x(k) = C*X;
    zmp_y(k) = C*Y;
    
    zmp_x_old = zmp_x(k);
    zmp_y_old = zmp_y(k);
    
    com_x(k) = X_new(0);
    com_y(k) = Y_new(0);
    
    u_sum_p_x = 0;
    u_sum_p_y = 0;
}

void setZmpref(int i){ // i = 0, 1, ... 13999

    mytime = i * 0.001;


    if(mytime < 3.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = 0.0;
    }
    else if(mytime < 4.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = 0.05;
    }
    else if(mytime < 5.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = -0.05;
    }
    else if(mytime < 6.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = 0.05;
    }
    else if(mytime < 7.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = -0.05;
    }
    else if(mytime < 8.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = 0.05;
    }
    else if(mytime < 9.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = -0.05;
    }
    else if(mytime < 10.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = 0.05;
    }
    else if(mytime < 11.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = -0.05;
    }
    else if(mytime < 12.0){
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = 0.05;
    }
    else{
    zmp_x_ref(i) = 0;
    zmp_y_ref(i) = 0;
    // zmp_y_ref(i) = 0;
    }
}

// for(int i = 0 ; i <= 15000 ; i++){
// mytime = i * dt;
 
// if(mytime < 3.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.0;
// }
// else if(mytime < 4.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.005;
// }
// else if(mytime < 5.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.01;
// }
// else if(mytime < 6.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.015;
// }
// else if(mytime < 7.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.02;
// }
// else if(mytime < 8.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.025;
// }
// else if(mytime < 9.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.03;
// }
// else if(mytime < 10.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.035;
// }
// else if(mytime < 11.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.04;
// }
// else if(mytime < 12.0){
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.045;
// }
// else{
// zmp_x_ref(i) = 0;
// zmp_y_ref(i) = 0.05;
// }
// }


typedef struct{
    double x;
    double y;
    double z;
} position;

position pos[3];

enum{
    BA, LF, RF
};

MatrixXd identity_matrix = MatrixXd::Identity(3,3);




namespace base{
    void setPosParam(double x=0, double y= 0, double z=0.25){
        pos[BA].x = x;
        pos[BA].y = y;
        pos[BA].z = z;
    }
}

namespace leftFoot{
    void setPosParam(double x=0, double y = 0.05, double z=0){
        pos[LF].x = x;
        pos[LF].y = y;
        pos[LF].z = z;
    }
}


namespace rightFoot{
    void setPosParam(double x=0, double y= -0.05, double z=0){  
        pos[RF].x = x;
        pos[RF].y = y;
        pos[RF].z = z;
    }
}





MatrixXd rotMatX(double q){
    MatrixXd tmp_m(3,3);
    tmp_m << 1, 0, 0,\
    0, cos(q), -sin(q),\
    0, sin(q), cos(q);
    return tmp_m;
}

MatrixXd rotMatY(double q){
    MatrixXd tmp_m(3,3);
    tmp_m << cos(q), 0, sin(q),\
    0 , 1, 0 ,\
    -sin(q), 0, cos(q);
    return tmp_m;
}

MatrixXd rotMatZ(double q){
    MatrixXd tmp_m(3,3);
    tmp_m << cos(q), -sin(q), 0,\
    sin(q), cos(q), 0, \
    0, 0, 1;
    return tmp_m;
}

int sign(double a){
    if(a>=0) return 1;
        return -1;
}

double func_1_cos(double t, double init, double final, double T)
{
    double des;
    
    des = init + ((final - init)/2) * (1-cos(PI*t/T));
    
    return des;
}






MatrixXd getTransformI0()
{
    // Frame 0 to Frame I

    MatrixXd tmp_m(4,4);

    // Alternative representations
    // 1. tmp_m = MatrixXd::Identity(4,4);
    // 2. tmp_m(m-1,n-1) = m,n element where m,n>=1;

    tmp_m << 1, 0, 0, 0, \
    0, 1, 0, 0, \
    0, 0, 1, 0, \
    0, 0, 0, 1;

    return tmp_m;
}

MatrixXd getTransform6E()
{
    // Frame E to Frame 6

    MatrixXd tmp_m(4,4);

    tmp_m << 1, 0, 0, 0, \
    0, 1, 0, 0, \
    0, 0, 1, -0.037, \
    0, 0, 0, 1;

    return tmp_m;
}

MatrixXd jointToTransform01(VectorXd q)
{
    // Frame 1 to Frame 0
    // q: generalized coordinates, q = [q1;q2;q3;q4;q5;q6]

    MatrixXd tmp_m(4,4);
    double qq = q(0);

    tmp_m << cos(qq), -sin(qq), 0, 0, \
    sin(qq), cos(qq), 0, 0.05, \
    0, 0, 1, 0, \
    0, 0, 0, 1;

    return tmp_m;
}

MatrixXd jointToTransform01_R(VectorXd q)
{
    // Frame 1 to Frame 0
    // q: generalized coordinates, q = [q1;q2;q3;q4;q5;q6]

    MatrixXd tmp_m(4,4);
    double qq = q(0);

    tmp_m << cos(qq), -sin(qq), 0, 0, \
    sin(qq), cos(qq), 0, -0.05, \
    0, 0, 1, 0, \
    0, 0, 0, 1;

    return tmp_m;
}

MatrixXd jointToTransform12(VectorXd q)
{
    // Frame 2 to Frame 1

    MatrixXd tmp_m(4,4);
    double qq = q(1);

    tmp_m << 1, 0, 0, 0, \
    0, cos(qq), -sin(qq), 0, \
    0, sin(qq), cos(qq), 0, \
    0 , 0, 0 , 1;

    return tmp_m;
}

MatrixXd jointToTransform23(VectorXd q)
{
    // Frame 3 to Frame 2
    MatrixXd tmp_m(4,4);
    double qq = q(2);

    tmp_m << cos(qq), 0, sin(qq), 0, \
    0 , 1, 0 , 0, \
    -sin(qq), 0, cos(qq), 0, \
    0 , 0, 0 , 1;

    return tmp_m;
}

MatrixXd jointToTransform34(VectorXd q)
{
    // Frame 4 to Frame 3
    MatrixXd tmp_m(4,4);
    double qq = q(3);

    tmp_m << cos(qq), 0, sin(qq), 0, \
    0 , 1, 0 , 0, \
    -sin(qq), 0, cos(qq), -0.133, \
    0 , 0, 0 , 1;

    return tmp_m;
}

MatrixXd jointToTransform45(VectorXd q)
{
    // Frame 5 to Frame 4
    MatrixXd tmp_m(4,4);
    double qq = q(4);

    tmp_m << cos(qq), 0, sin(qq), 0, \
    0 , 1, 0 , 0, \
    -sin(qq), 0, cos(qq), -0.138, \
    0 , 0, 0 , 1;

    return tmp_m;
}

MatrixXd jointToTransform56(VectorXd q)
{
    // Frame 5 to Frame 4
    MatrixXd tmp_m(4,4);
    double qq = q(5);

    tmp_m << 1, 0, 0, 0, \
    0, cos(qq), -sin(qq), 0, \
    0, sin(qq), cos(qq), 0, \
    0 , 0, 0 , 1;

    return tmp_m;
}

MatrixXd jointToPosition(VectorXd q, MatrixXd GB)
{
    // Extract position vector(3x1)

    VectorXd tmp_v = VectorXd::Zero(3);
    MatrixXd tmp_m(4,4);

    tmp_m = GB*
    jointToTransform01(q)*
    jointToTransform12(q)*
    jointToTransform23(q)*
    jointToTransform34(q)*
    jointToTransform45(q)*
    jointToTransform56(q)*
    getTransform6E();

    // tmp_v = tmp_m.block(0,3,3,1);

    tmp_v(0) = tmp_m(0,3);
    tmp_v(1) = tmp_m(1,3);
    tmp_v(2) = tmp_m(2,3);

    return tmp_v;
}

MatrixXd jointToPosition_R(VectorXd q, MatrixXd GB)
{
    // Extract position vector(3x1)

    VectorXd tmp_v = VectorXd::Zero(3);
    MatrixXd tmp_m(4,4);

    tmp_m = GB*
    jointToTransform01_R(q)*
    jointToTransform12(q)*
    jointToTransform23(q)*
    jointToTransform34(q)*
    jointToTransform45(q)*
    jointToTransform56(q)*
    getTransform6E();

    // tmp_v = tmp_m.block(0,3,3,1);

    tmp_v(0) = tmp_m(0,3);
    tmp_v(1) = tmp_m(1,3);
    tmp_v(2) = tmp_m(2,3);

    return tmp_v;
}




VectorXd Geometric_IK_L(){
    double A, B, C; // thigh, calf
    double alpha, c5;
    double tmp, theta_prime;
    double r_x, r_y, r_z;
    VectorXd B_r_BH(3) ;
    VectorXd r_GB(3), r_GH(3), r_GA(3), r_GF(3), r_FA(3), r(3);
    VectorXd q(6);
    MatrixXd C_GB(3,3), C_GH(3,3), C_GA(3,3), C_GF(3,3);
    
    B_r_BH << 0, 0.05, 0;
    
    // C_GB = T_GB.block(0,0,3,3);
    // C_GF = T_GF.block(0,0,3,3);
    C_GB = identity_matrix;
    C_GF = identity_matrix;
    // r_GF = T_GF.block(0,3,3,1);
    r_GF << pos[LF].x, pos[LF].y, pos[LF].z;
    r_FA << 0, 0, 0.037 ;

    r_GA = r_GF + r_FA ;
    
    A = 0.133;
    B = 0.138;
    
    // r_GB << T_GB(0,3), T_GB(1,3), T_GB(2,3);
    r_GB << pos[BA].x, pos[BA].y, pos[BA].z;
    r_GH = r_GB + C_GB * B_r_BH;
    
    C_GA = C_GF;
    r= C_GA.transpose() * (r_GH - r_GA); // r = A_r_AH = C_AG * G_r_AH
    r_x = r(0);
    r_y = r(1);
    r_z = r(2);
    C= r.norm(); // r의 길이

    // c5 = cos((A*A + B*B - C*C) / (2 * A * B));
    
    // if(c5 >= 1) q(3) = 0.0;
    
    // else if(c5 <= -1) q(3) = PI;
    
    // else q(3)=-acos((A*A + B*B - C*C) / (2 * A * B)) + PI;
    
    tmp = (A*A + B*B - C*C) / (2 * A * B);
    // acos함수의 정의역은 -1부터 1 사이어야 한다.
    if(tmp >= 1)
    theta_prime = 0.0;
    else if(tmp <= -1)
    theta_prime = PI;
    else
    theta_prime = acos(tmp);

    q(3)=PI - theta_prime; // knee pitch

    alpha = asin((A*sin(PI - q(3)))/C);
    
    q(5) = atan2(r_y, r_z); // ankle roll
    q(4) = -atan2(r_x, sign(r_z)*sqrt(r_y*r_y + r_z*r_z)) - alpha; // ankle pitch
    
    C_GH = C_GB.transpose() * C_GF * rotMatX(-q(5)) * rotMatY(-q(3)-q(4));
    
    q(0) = atan2(-C_GH(0,1), C_GH(1,1)); // hip yaw
    q(1) = atan2(C_GH(2,1), -sin(q(0)) * C_GH(0,1) + C_GH(1,1) * cos(q(0))); // hip roll
    q(2) = atan2(-C_GH(2,0), C_GH(2,2)); // hip pitch
    
    return q;
 
}


VectorXd Geometric_IK_R(){
    double A, B, C; // thigh, calf
    double alpha, c5;
    double tmp, theta_prime;
    double r_x, r_y, r_z;
    VectorXd B_r_BH(3) ;
    VectorXd r_GB(3), r_GH(3), r_GA(3), r_GF(3), r_FA(3), r(3);
    VectorXd q(6);
    MatrixXd C_GB(3,3), C_GH(3,3), C_GA(3,3), C_GF(3,3);
    
    B_r_BH << 0, -0.05, 0; // 왼쪽과 오른쪽이 다른 부분이다.
    
    // C_GB = T_GB.block(0,0,3,3);
    // C_GF = T_GF.block(0,0,3,3);
    C_GB = identity_matrix;
    C_GF = identity_matrix;
    // r_GF = T_GF.block(0,3,3,1);
    r_GF << pos[RF].x, pos[RF].y, pos[RF].z;
    r_FA << 0, 0, 0.037 ;
    r_GA = r_GF + r_FA ;
    
    A = 0.133;
    B = 0.138;
    
    // r_GB << T_GB(0,3), T_GB(1,3), T_GB(2,3);
    r_GB << pos[BA].x, pos[BA].y, pos[BA].z;
    r_GH = r_GB + C_GB * B_r_BH;
    
    C_GA = C_GF;
    r= C_GA.transpose() * (r_GH - r_GA); // r = A_r_AH = C_AG * G_r_AH
    r_x = r(0);
    r_y = r(1);
    r_z = r(2);
    C= r.norm(); // r의 길이

    // c5 = cos((A*A + B*B - C*C) / (2 * A * B));
    
    // if(c5 >= 1) q(3) = 0.0;
    
    // else if(c5 <= -1) q(3) = PI;
    
    // else q(3)=-acos((A*A + B*B - C*C) / (2 * A * B)) + PI;
    
    tmp = (A*A + B*B - C*C) / (2 * A * B);
    // acos함수의 정의역은 -1부터 1 사이어야 한다.
    if(tmp >= 1)
    theta_prime = 0.0;
    else if(tmp <= -1)
    theta_prime = PI;
    else
    theta_prime = acos(tmp);

    q(3)=PI - theta_prime; // knee pitch

    alpha = asin((A*sin(PI - q(3)))/C);
    
    q(5) = atan2(r_y, r_z); // ankle roll
    q(4) = -atan2(r_x, sign(r_z)*sqrt(r_y*r_y + r_z*r_z)) - alpha; // ankle pitch
    
    C_GH = C_GB.transpose() * C_GF * rotMatX(-q(5)) * rotMatY(-q(3)-q(4));
    
    q(0) = atan2(-C_GH(0,1), C_GH(1,1)); // hip yaw
    q(1) = atan2(C_GH(2,1), -sin(q(0)) * C_GH(0,1) + C_GH(1,1) * cos(q(0))); // hip roll
    q(2) = atan2(-C_GH(2,0), C_GH(2,2)); // hip pitch
    
    return q;
    
}






void GetJoint(const std_msgs::Float64MultiArray &msg);

ros::Subscriber GET_JOINT;


namespace gazebo{
    class PIDJoints : public ModelPlugin
    {
        double dt;

        double LP_ANGLE = 0;
        double LPm_ANGLE = 0;
        double LPd_ANGLE = 0;
        double LK_ANGLE = 0;
        double LA_ANGLE = 0;
        double LF_ANGLE = 0;
        
        double RP_ANGLE = 0;
        double RPm_ANGLE = 0;
        double RPd_ANGLE = 0;
        double RK_ANGLE = 0;
        double RA_ANGLE = 0;
        double RF_ANGLE = 0;
        

        double LP_ANGLE_er = 0;
        double LPm_ANGLE_er = 0;
        double LPd_ANGLE_er = 0;
        double LK_ANGLE_er = 0;
        double LA_ANGLE_er = 0;
        double LF_ANGLE_er = 0;
        
        double RP_ANGLE_er = 0;
        double RPm_ANGLE_er = 0;
        double RPd_ANGLE_er = 0;
        double RK_ANGLE_er = 0;
        double RA_ANGLE_er = 0;
        double RF_ANGLE_er = 0;

        physics::LinkPtr base_link ;
        physics::LinkPtr L_P_link ;
        physics::LinkPtr L_Pm_link ;
        physics::LinkPtr L_Pd_link ;
        physics::LinkPtr L_K_link ;
        physics::LinkPtr L_A_link ;
        physics::LinkPtr L_F_link ;
        
        physics::LinkPtr R_P_link ;
        physics::LinkPtr R_Pm_link ;
        physics::LinkPtr R_Pd_link ;
        physics::LinkPtr R_K_link ;
        physics::LinkPtr R_A_link ;
        physics::LinkPtr R_F_link ;
        
        
        physics::JointPtr LP ;
        physics::JointPtr LPm ;
        physics::JointPtr LPd ;
        physics::JointPtr LK ;
        physics::JointPtr LA ;
        physics::JointPtr LF ;
        
        
        physics::JointPtr RP ;
        physics::JointPtr RPm ;
        physics::JointPtr RPd ;
        physics::JointPtr RK ;
        physics::JointPtr RA ;
        physics::JointPtr RF ;
        

        
        common::PID LPpid ;
        common::PID LPmpid ;
        common::PID LPdpid ;
        common::PID LKpid ;
        common::PID LApid ;
        common::PID LFpid ;
        
        common::PID RPpid ;
        common::PID RPmpid ;
        common::PID RPdpid ;
        common::PID RKpid ;
        common::PID RApid ;
        common::PID RFpid ;
        
        
        physics::ModelPtr model;

        common::Time last_update_time;
        event::ConnectionPtr update_connection_;

        ros::NodeHandle nh;
        ros::Publisher P_Times;


        ros::Publisher zmp_x_ref_pub;
        ros::Publisher zmp_y_ref_pub;
        ros::Publisher com_x_pub;
        ros::Publisher com_y_pub;

        std_msgs::Float64 zmp_x_ref_msg;
        std_msgs::Float64 zmp_y_ref_msg;
        std_msgs::Float64 com_x_msg;
        std_msgs::Float64 com_y_msg;

        
        public:
        void Load(physics::ModelPtr _model, sdf::ElementPtr /*_sdf*/);
        void UpdatePID();
        int rqt(int argc, char **argv);

        void publishtopic(int i);
    
    };

GZ_REGISTER_MODEL_PLUGIN(PIDJoints);
 
}

void gazebo::PIDJoints::Load(physics::ModelPtr _model, sdf::ElementPtr /*_sdf*/)
{
    cout << "Load Function" << endl;

    getGi();
    getGx();
    getGp();

    //////////////// 파일 출력

    fout.open("/home/ubuntu/catkin_ws/src/kubot22_pkgs/src/com_y_2.txt");
    ROS_INFO("open fout");

    ///////////////// 토픽 퍼블리셔

    zmp_x_ref_pub = nh.advertise<std_msgs::Float64>("zmp_x_ref", 1000);
    zmp_y_ref_pub = nh.advertise<std_msgs::Float64>("zmp_y_ref", 1000);

    com_x_pub = nh.advertise<std_msgs::Float64>("com_x", 1000);
    com_y_pub = nh.advertise<std_msgs::Float64>("com_y", 1000);


    ////////////// 정의

    i=0;

    All_time_trajectory = 15.0;
    mydt = 0.001;
    // N = 300;

    n = 15001;

    z_c = 0.25;
    g = 9.81;

    A << 1, mydt, mydt*mydt/2, \
    0, 1, mydt, \
    0, 0, 1 ;

    B << mydt*mydt*mydt/6, mydt*mydt/2, mydt;

    C << 1, 0, -z_c/g;

    // Gi =    8.8764998e+02;

    Gx <<   Gx_1, Gx_2, Gx_3;

    ////// Gp 즉, 배열은 선언과 정의를 동시에 해야한다.

    X = VectorXd::Zero(3); //3*1

    X_new = VectorXd::Zero(3);

    zmp_x_ref = RowVectorXd::Zero(n);

    u_x = 0;

    Y = VectorXd::Zero(3);

    Y_new = VectorXd::Zero(3);

    zmp_y_ref = RowVectorXd::Zero(n);

    u_y = 0;



    zmp_x = RowVectorXd::Zero(n);

    zmp_x_old = 0;

    zmp_y = RowVectorXd::Zero(n);

    zmp_y_old = 0;

    com_x = RowVectorXd::Zero(n);

    com_y = RowVectorXd::Zero(n);

    err_x = 0;
    err_y = 0;

    sum_e_x = 0;
    sum_e_y = 0;

    u_sum_p_x = 0;
    u_sum_p_y = 0;

    ///////////////////

    // for(int i =0 ; i <= 15000 ; i+= 1000){
    // ROS_INFO("com_x(%d) : %f", i, com_x(i));
    // }

    // for(int i =0 ; i <= 15000 ; i+= 1000){
    // ROS_INFO("com_y(%d) : %f", i, com_y(i));
    // }

    //////////////////

    this->model = _model;
    
    model = _model;
    
    this->base_link = this->model->GetLink("base_link");
    
    this->L_P_link = this->model->GetLink("L_P_link");
    this->L_Pm_link = this->model->GetLink("L_Pm_link");
    this->L_Pd_link = this->model->GetLink("L_Pd_link");
    this->L_K_link = this->model->GetLink("L_K_link");
    this->L_A_link = this->model->GetLink("L_A_link");
    this->L_F_link = this->model->GetLink("L_F_link");
    
    this->R_P_link = this->model->GetLink("R_P_link");
    this->R_Pm_link = this->model->GetLink("R_Pm_link");
    this->R_Pd_link = this->model->GetLink("R_Pd_link");
    this->R_K_link = this->model->GetLink("R_K_link");
    this->R_A_link = this->model->GetLink("R_A_link");
    this->R_F_link = this->model->GetLink("R_F_link");
    
    
    this->LP = this->model->GetJoint("LP");
    this->LPm = this->model->GetJoint("LPm");
    this->LPd = this->model->GetJoint("LPd");
    this->LK = this->model->GetJoint("LK");
    this->LA = this->model->GetJoint("LA");
    this->LF = this->model->GetJoint("LF");
    
    this->RP = this->model->GetJoint("RP");
    this->RPm = this->model->GetJoint("RPm");
    this->RPd = this->model->GetJoint("RPd");
    this->RK = this->model->GetJoint("RK");
    this->RA = this->model->GetJoint("RA");
    this->RF = this->model->GetJoint("RF");
    

    
    #if GAZEBO_MAJOR_VERSION >= 8
    last_update_time = model->GetWorld()->SimTime();
    #else
    last_update_time = model->GetWorld()->GetSimTime();
    #endif
    this->update_connection_ = event::Events::ConnectWorldUpdateBegin(boost::bind(&PIDJoints::UpdatePID, this));
    
    // p gain d gain
    this->LPpid.Init( 80, 0, 0.20, 200, -200, 1000, -1000);
    this->LPmpid.Init( 90, 0, 0.30, 200, -200, 1000, -1000);
    this->LPdpid.Init(110, 0, 0.80, 200, -200, 1000, -1000);
    this->LKpid.Init( 110, 0, 0.70, 200, -200, 1000, -1000);
    this->LApid.Init( 110, 0, 0.80, 200, -200, 1000, -1000);
    this->LFpid.Init( 90, 0, 0.30, 200, -200, 1000, -1000);
    
    this->RPpid.Init( 80, 0, 0.20, 200, -200, 1000, -1000);
    this->RPmpid.Init( 90, 0, 0.30, 200, -200, 1000, -1000);
    this->RPdpid.Init(110, 0, 0.80, 200, -200, 1000, -1000);
    this->RKpid.Init( 110, 0, 0.70, 200, -200, 1000, -1000);
    this->RApid.Init( 110, 0, 0.80, 200, -200, 1000, -1000);
    this->RFpid.Init( 90, 0, 0.30, 200, -200, 1000, -1000);
    
    
    
    int version_test;
    
    version_test = rbdl_get_api_version();
    
    printf("RBDL API version = %d\n", version_test);
    
    Addons::URDFReadFromFile("/home/ubuntu/.gazebo/models/kubot22/urdf/kubot22.urdf", kubot22_Model, false, false);
    
    kubot22_Model -> gravity = Eigen::Vector3d(0., 0., -9.81);
    
    RobotState = VectorNd::Zero(kubot22_Model->q_size);
    RobotStatedot = VectorNd::Zero(kubot22_Model->q_size);
    RobotStatedot2 = VectorNd::Zero(kubot22_Model->q_size);

    
    GET_JOINT = nh.subscribe("desQ_array", 100, GetJoint);
    
    
    ros::Rate loop_rate(1000);
    
    
    q_init << 0, 0, 0, 0, 0, 0;
    q_init_R << 0, 0, 0, 0, 0, 0;





 
}

void gazebo::PIDJoints::UpdatePID()//여러번 실행
{
    //cout << "A" << endl;
    #if GAZEBO_MAJOR_VERSION >= 8
    common::Time current_time = model->GetWorld()->SimTime();
    #else
    common::Time current_time = model->GetWorld()->GetSimTime(); // get simulation time
    #endif
    dt = current_time.Double() - this->last_update_time.Double();

    cnt = cnt + dt;
    //cnt=0;

    
    this->last_update_time = current_time;

    // 현재의 각도
    LP_ANGLE = this-> LP ->Position(2);
    LPm_ANGLE = this-> LPm ->Position(1);
    LPd_ANGLE = this-> LPd ->Position(1);
    LK_ANGLE = this-> LK ->Position(1);
    LA_ANGLE = this-> LA ->Position(1);
    LF_ANGLE = this-> LF ->Position(1);
    
    RP_ANGLE = this-> RP ->Position(2);
    RPm_ANGLE = this-> RPm ->Position(1);
    RPd_ANGLE = this-> RPd ->Position(1);
    RK_ANGLE = this-> RK ->Position(1);
    RA_ANGLE = this-> RA ->Position(1);
    RF_ANGLE = this-> RF ->Position(1);
    
    // 아무동작도 안하고 베이스 z만 0.25로 움직였을 때 ROS_INFO("LPd_ANGLE : %f", LPd_ANGLE);
    // 로 출력한 값들
    q_cal << 0, 0, -0.68, 1.3, -0.65, 0;
    q_cal_R << 0, 0, -0.68, 1.3, -0.65, 0;


    if (cnt < 3 ) {
        ROS_INFO("port status : %d", fout.is_open());
        // ROS_INFO("LPd_ANGLE : %f", LPd_ANGLE);
        // ROS_INFO("LK_ANGLE : %f", LK_ANGLE);
        // ROS_INFO("LA_ANGLE : %f", LA_ANGLE);

        LP_ANGLE_er = LP_ANGLE - func_1_cos(cnt, q_init(0), q_cal(0), T);
        LPm_ANGLE_er = LPm_ANGLE - func_1_cos(cnt, q_init(1), q_cal(1), T);
        LPd_ANGLE_er = LPd_ANGLE - func_1_cos(cnt, q_init(2), q_cal(2), T);
        LK_ANGLE_er = LK_ANGLE - func_1_cos(cnt, q_init(3), q_cal(3), T);
        LA_ANGLE_er = LA_ANGLE - func_1_cos(cnt, q_init(4), q_cal(4), T);
        LF_ANGLE_er = LF_ANGLE - func_1_cos(cnt, q_init(5), q_cal(5), T);
        
        RP_ANGLE_er = RP_ANGLE - func_1_cos(cnt, q_init_R(0), q_cal_R(0), T);
        RPm_ANGLE_er = RPm_ANGLE - func_1_cos(cnt, q_init_R(1), q_cal_R(1), T);
        RPd_ANGLE_er = RPd_ANGLE - func_1_cos(cnt, q_init_R(2), q_cal_R(2), T);
        RK_ANGLE_er = RK_ANGLE - func_1_cos(cnt, q_init_R(3), q_cal_R(3), T);
        RA_ANGLE_er = RA_ANGLE - func_1_cos(cnt, q_init_R(4), q_cal_R(4), T);
        RF_ANGLE_er = RF_ANGLE - func_1_cos(cnt, q_init_R(5), q_cal_R(5), T);
    }

    else if (cnt >= 3 && cnt < 17) {
        cout << "i = " << i << endl;

        setZmpref(i);
        cout << "zmp_y_ref = " << zmp_y_ref(i) << endl;

        // fout << zmp_y_ref(i) << endl;

        setCom(i);
        cout << "com_y = " << com_y(i) << endl;

        fout << com_y(i) << endl;
        // fout << zmp_y(i) << endl;


        publishtopic(i);

        base::setPosParam(0, com_y(i), 0.25);
        rightFoot::setPosParam();
        leftFoot::setPosParam();
        
        q_R = Geometric_IK_R();
        q = Geometric_IK_L();

        LP_ANGLE_er = LP_ANGLE - q(0);
        LPm_ANGLE_er = LPm_ANGLE - q(1);
        LPd_ANGLE_er = LPd_ANGLE - q(2);
        LK_ANGLE_er = LK_ANGLE - q(3);
        LA_ANGLE_er = LA_ANGLE - q(4);
        LF_ANGLE_er = LF_ANGLE - q(5);

        RP_ANGLE_er = RP_ANGLE - q_R(0);
        RPm_ANGLE_er = RPm_ANGLE - q_R(1);
        RPd_ANGLE_er = RPd_ANGLE - q_R(2);
        RK_ANGLE_er = RK_ANGLE - q_R(3);
        RA_ANGLE_er = RA_ANGLE - q_R(4);
        RF_ANGLE_er = RF_ANGLE - q_R(5);
        
        i++;
        if( i >= 14702 )
            cnt = 17.001;
    }

    else if (cnt >= 17 && cnt < 20) {
        ROS_INFO("cnt : %f", cnt);
        base::setPosParam(0, 0, 0.25);
        rightFoot::setPosParam();
        // leftFoot::setPosParam(0, 0, func_1_cos(cnt-17, 0, 0.01, 3));
        leftFoot::setPosParam();
        
        q_R = Geometric_IK_R();
        q = Geometric_IK_L();

        LP_ANGLE_er = LP_ANGLE - q(0);
        LPm_ANGLE_er = LPm_ANGLE - q(1);
        LPd_ANGLE_er = LPd_ANGLE - q(2);
        LK_ANGLE_er = LK_ANGLE - q(3);
        LA_ANGLE_er = LA_ANGLE - q(4);
        LF_ANGLE_er = LF_ANGLE - q(5);

        RP_ANGLE_er = RP_ANGLE - q_R(0);
        RPm_ANGLE_er = RPm_ANGLE - q_R(1);
        RPd_ANGLE_er = RPd_ANGLE - q_R(2);
        RK_ANGLE_er = RK_ANGLE - q_R(3);
        RA_ANGLE_er = RA_ANGLE - q_R(4);
        RF_ANGLE_er = RF_ANGLE - q_R(5);

    }

    else{
        if(fout.is_open() == true){
            ROS_INFO("close fout ");
            fout.close();
        }
        //여기 부분이 끝나는 각도와 맞아야 한다.
        LP_ANGLE_er = LP_ANGLE - q(0);
        LPm_ANGLE_er = LPm_ANGLE - q(1);
        LPd_ANGLE_er = LPd_ANGLE - q(2);
        LK_ANGLE_er = LK_ANGLE - q(3);
        LA_ANGLE_er = LA_ANGLE - q(4);
        LF_ANGLE_er = LF_ANGLE - q(5);
        
        RP_ANGLE_er = RP_ANGLE - q_R(0);
        RPm_ANGLE_er = RPm_ANGLE - q_R(1);
        RPd_ANGLE_er = RPd_ANGLE - q_R(2);
        RK_ANGLE_er = RK_ANGLE - q_R(3);
        RA_ANGLE_er = RA_ANGLE - q_R(4);
        RF_ANGLE_er = RF_ANGLE - q_R(5);
    }

    this -> LPpid.Update(LP_ANGLE_er, dt);
    this -> LPmpid.Update(LPm_ANGLE_er, dt);
    this -> LPdpid.Update(LPd_ANGLE_er, dt);
    this -> LKpid.Update(LK_ANGLE_er, dt);
    this -> LApid.Update(LA_ANGLE_er, dt);
    this -> LFpid.Update(LF_ANGLE_er, dt);
    
    this -> RPpid.Update(RP_ANGLE_er, dt);
    this -> RPmpid.Update(RPm_ANGLE_er, dt);
    this -> RPdpid.Update(RPd_ANGLE_er, dt);
    this -> RKpid.Update(RK_ANGLE_er, dt);
    this -> RApid.Update(RA_ANGLE_er, dt);
    this -> RFpid.Update(RF_ANGLE_er, dt);


    
    for (int LnJoint = 0; LnJoint < 12; LnJoint++) {
        RobotState(LnJoint) = NowDq(LnJoint);
        RobotStatedot(LnJoint) = NowDq(LnJoint);
        RobotStatedot2(LnJoint) = NowD2q(LnJoint);
    }

    
    LP ->SetForce(2, LPpid.GetCmd()); //setForce(axis,Force value)
    LPm->SetForce(1, LPmpid.GetCmd());
    LPd->SetForce(1, LPdpid.GetCmd());
    LK ->SetForce(1, LKpid.GetCmd());
    LA ->SetForce(1, LApid.GetCmd());
    LF ->SetForce(1, LFpid.GetCmd());

    RP ->SetForce(2, RPpid.GetCmd()); //setForce(axis,Force value)
    RPm->SetForce(1, RPmpid.GetCmd());
    RPd->SetForce(1, RPdpid.GetCmd());
    RK ->SetForce(1, RKpid.GetCmd());
    RA ->SetForce(1, RApid.GetCmd());
    RF ->SetForce(1, RFpid.GetCmd());
}

void GetJoint(const std_msgs::Float64MultiArray &msg)
{
    LJ[0] = msg.data[0];
    LJ[1] = msg.data[1];
    LJ[2] = msg.data[2];
    LJ[3] = msg.data[3];
    LJ[4] = msg.data[4];
    LJ[5] = msg.data[5];

    RJ[0] = msg.data[6];
    RJ[1] = msg.data[7];
    RJ[2] = msg.data[8];
    RJ[3] = msg.data[9];
    RJ[4] = msg.data[10];
    RJ[5] = msg.data[11];
}


void gazebo::PIDJoints::publishtopic(int i){
    zmp_x_ref_msg.data = zmp_x(i);
    zmp_y_ref_msg.data = zmp_y(i);
    com_x_msg.data = com_x(i);
    com_y_msg.data = com_y(i);

    zmp_x_ref_pub.publish(zmp_x_ref_msg);
    zmp_y_ref_pub.publish(zmp_y_ref_msg);
    com_x_pub.publish(com_x_msg);
    com_y_pub.publish(com_y_msg);
}
