#include <iostream>
#include <math.h>
using namespace std;

#define h0       1                // начальный шаг
#define alpha0   1                // начальный угол delta_alpha
#define v        10               // скорость АЭ волны в ОК

/* координаты датчиков */
#define x1       -1
#define y1       -2
#define x2       3
#define y2       -1
#define x3       0.5
#define y3       3
/*--------------------*/

#define min_f1   0.1             // погрешность нахождения невязки f1
#define min_f2   0.1             // погрешность нахождения невязки f2
#define max_i    50               // максимальное число итераций сходимости невязки f1
#define max_j    200              // максимальное число итераций сходимости невязки f2



double t21, t32, t31, Xx12, Yx12, angle_beta, Xc, Yc, buf_x, buf_y, a, b, S12, i, j, S23, S13, x12, y12;
double h[3], alpha[2], f1[2], f2[2], delta_alpha[3], f1_b[2];

int main()
{
    while (1)
    {

    S12=sqrt(pow((x2-x1),2)+pow((y2-y1),2));
    S23=sqrt(pow((x3-x2),2)+pow((y3-y2),2));
    S13=sqrt(pow((x3-x1),2)+pow((y3-y1),2));

    /* Ввод значений разности времен прихода сигналов АЭ на датчики */
    cout << "Enter t21, t32, t31: ";
    cin >> t21 >> t32 >> t31;
    /*-------------------------------------------------------------*/

    /* Вычисление 0-й невязки f1 и f2 */
    angle_beta=M_PI-atan(y1/x1)-acos((sqrt(pow(x1,2)+pow(y1,2))+pow((x2-x1),2)+pow((y2-y1),2)-sqrt(pow(x2,2)+pow(y2,2)))/(2*sqrt(pow(x1,2)+pow(y1,2))*S12));
    Xx12=(S12-t21*v)/2*cos(angle_beta)+x1;                                  // абсцисса вершины гиперболы
    Yx12=((sqrt(pow((x2-x1),2)+pow((y2-y1),2))-t21*v)/2)*sin(angle_beta)+y1; // ордината вершины гиперболы
    alpha[0]=M_PI-angle_beta;                                              // угол наклона линии, соединяющей два датчика, с осью абсцисс
    f1[0]=(S12-t21*v)/2-t21*v;                                               // 0-ая невязка f1
    f2[0]=sqrt(pow((x3-Xx12),2)+pow((y3-Yx12),2))-(S12-t21*v)/2-t31*v;       // 0-ая невязка f2
    Xc=S12*cos(alpha[0])+x1;                                               // абсцисса центра гиперболы
    Yc=S12*sin(alpha[0])+y1;                                               // ордината центра гиперболы
    /*-------------------------------*/

    /* откладываем единичный отрезок и проверяем правильность выбора направления */
    buf_x=(Xx12-Xc)*cos(alpha[0])+(Yx12-Yc)*sin(alpha[0]);
    buf_y=(Yx12-Xc)*cos(alpha[0])-(Xx12-Yc)*sin(alpha[0]);
    alpha[1]=alpha[0]+atan(buf_x/buf_y*pow(S12,2)/pow((v*t21),2)-1); // угол между касательной к вершине гиперболы и осью абсцисс
    buf_x=h0*cos(alpha[1])+Xx12;
    buf_y=h0*sin(alpha[1])+Yx12;
    f1[1]=sqrt(pow((x2-buf_x),2)+pow((y2-buf_y),2))-sqrt(pow((x1-buf_x),2)+pow((y1-buf_y),2))-t21*v;
    f2[1]=sqrt(pow((x3-buf_x),2)+pow((y3-buf_y),2))-sqrt(pow((x1-buf_x),2)+pow((y1-buf_y),2))-t31*v;
    if (fabs(f2[1])>fabs(f2[0]))                                                       // если направление выбрано неверно
    {
        alpha[1]=alpha[1]-M_PI;                                             // меняем направление единичного отрезка
        buf_x=h0*cos(alpha[1])+Xx12;
        buf_y=h0*sin(alpha[1])+Yx12;
        f1[1]=sqrt(pow((x2-buf_x),2)+pow((y2-buf_y),2))-sqrt(pow((x1-buf_x),2)+pow((y1-buf_y),2))-t21*v; // 1-я невязка f1
        f2[1]=sqrt(pow((x3-buf_x),2)+pow((y3-buf_y),2))-sqrt(pow((x1-buf_x),2)+pow((y1-buf_y),2))-t31*v; // 1-я невязка f2
    }
    /*---------------------------------------------------------------------------*/

    h[0]=0;
    h[1]=h0;
    cout<<"f1[1]="<<f1[1]<<"  f2[1]="<<f2[1]<<endl;
    while (((fabs(f1[1])>min_f1)||(fabs(f2[1])>min_f2))&&(j<=max_j))
    {
        delta_alpha[0]=alpha0;
        delta_alpha[1]=delta_alpha[0]/2;
        delta_alpha[2]=delta_alpha[0]/4;
        i=0;
        f1_b[0]=f1[0];
        f1_b[1]=f1[1];
        /* нахождение невязки f1 методом секущих */
        while ((fabs(f1_b[1])>min_f1)&&(i<=max_i))
        {
            a=(f1_b[1]-f1_b[0])/(delta_alpha[1]-delta_alpha[0]);
            b=(f1_b[0]*delta_alpha[1]-f1_b[1]*delta_alpha[0])/(delta_alpha[1]-delta_alpha[0]);
            f1_b[0]=f1_b[1];
            f1_b[1]=a*delta_alpha[2]+b;

            if (((f1_b[0]>0)&&(f1_b[1]<0))||((f1_b[0]<0)&&(f1_b[1]>0)))
            {
                delta_alpha[2]=delta_alpha[2]/2;
                f1_b[1]=a*delta_alpha[2]+b;
            }
            delta_alpha[0]=delta_alpha[1];
            delta_alpha[1]=delta_alpha[2];
            delta_alpha[2]=-(f1_b[0]*delta_alpha[1]-f1_b[1]*delta_alpha[0])/(f1_b[1]-f1_b[0]);
            i++;
        }
        /*--------------------------------------*/
        cout<<"i="<<i<<endl;
        if ((i>max_i)&&(fabs(f1[1])>min_f1))
        {
            break;
        }

        /* нахождение шага hi+1 */
        h[2]=-(f2[0]*h[1]-f2[1]*h[0])/(f2[1]-f2[0]);
        x12=(buf_x-Xc)*cos(alpha[0])+(buf_y-Yc)*sin(alpha[0]);
        y12=(buf_y-Xc)*cos(alpha[0])-(buf_x-Yc)*sin(alpha[0]);
        alpha[1]=atan(x12/y12*pow(S12,2)/pow((v*t21),2)-1)+alpha[0];
        /*----------------------*/

        /* вычисление текущей точки на гиперболе */
        buf_x=h[2]*cos(alpha[1]+delta_alpha[2]*M_PI/180)+buf_x;
        buf_y=h[2]*sin(alpha[1]+delta_alpha[2]*M_PI/180)+buf_y;
        /*--------------------------------------*/

        f1[0]=f1[1];
        f2[0]=f2[1];
        f1[1]=sqrt(pow((x2-buf_x),2)+pow((y2-buf_y),2))-sqrt(pow((x1-buf_x),2)+pow((y1-buf_y),2))-t21*v;
        f2[1]=sqrt(pow((x3-buf_x),2)+pow((y3-buf_y),2))-sqrt(pow((x1-buf_x),2)+pow((y1-buf_y),2))-t31*v;
        cout<<"f1[1]="<<f1[1]<<"  f2[1]="<<f2[1]<<endl;
        h[0]=h[1]; h[1]=h[2];
        j++;

    }
    cout<<"j="<<j<<endl;
    if (((j<=max_j+1)&&(fabs(f2[1])<=min_f2))&&(fabs(f1[1])<=min_f1))
    {
        cout << buf_x << ";" << buf_y<<endl;
    }
    j=0;
    }
    return 0;
}
