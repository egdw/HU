#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>


using namespace cv;
using namespace std;
//bool HuMoment(Mat* img);
double * getHu(Moments m);
vector<double> HuMoment(IplImage* img);
void calcCiCS(vector<vector<double *>> * datas,IplImage *channel_h,IplImage *channel_s,IplImage *channel_v);
int main() {
//    std::cout << "Hello, World!" << std::endl;
//    //efault is bgr
//    Mat input = imread("/home/egdw/图片/数据集/0001-01.png");
//    //so,we need to convert bgr to hsv and lab
////    Mat hsv;
//    Mat lab;
//    Mat gray;
////    cvtColor(input,hsv,COLOR_BGR2HSV);
//    cvtColor(input,lab,COLOR_BGR2Lab);
//    cvtColor(input,gray,COLOR_BGR2GRAY);
//
////    cout<< hsv.channels() << " " << hsv.at<Vec<float,3>>(2,2) << endl;
//    cout<< lab.channels() << " " << lab.at<Vec<float,3>>(0,0) << endl;
//
//    cout<< lab.rows  << lab.cols<< endl;
    //calculate hu
//    namedWindow("Test");
//    imshow("display",gray);
//    waitKey();

    IplImage * input = cvLoadImage("/home/egdw/图片/数据集/0001-01.png", 1);
    IplImage *hsv = cvCreateImage(cvSize(input->width, input->height), input->depth, input->nChannels);  //注意图像必须和输入图像的size，颜色位深度，通道一致
    cvZero(hsv);
    cvCvtColor(input, hsv, CV_BGR2HSV);
    IplImage *channel_h=cvCreateImage(cvSize(hsv->width,hsv->height),IPL_DEPTH_8U,1);
    IplImage *channel_s=cvCreateImage(cvSize(hsv->width,hsv->height),IPL_DEPTH_8U,1);
    IplImage *channel_v=cvCreateImage(cvSize(hsv->width,hsv->height),IPL_DEPTH_8U,1);
    cvSplit(input, channel_h, channel_s, channel_v, NULL);

    vector<vector<double *>> datas;
//    cout << "exe" << endl;

    calcCiCS(&datas,channel_h,channel_s,channel_v);
//    cout << "exe" << endl;
//    cout << datas.size() << endl;
//    for(int i = 0;i<datas.size();i++){
//        for(int j = 0;j<datas[i].size();j++){
//            cout << "(" << datas[i][j][0] << "," << datas[i][j][1] << ") ";
//        }
//        cout << endl;
//    }

//    calcCiCS();
//    calcCiCS(channel_v);


    //then
   //et the hu from h s v channel.
//    HuMoment(channel_h);
//    HuMoment(channel_s);
//    HuMoment(channel_v);
//
//    // to cal cicsH cicsS cicsV
//    cout << float(((uchar*)channel_s->imageData)[0*channel_s->widthStep+4]) << endl;
//    int i,j;
//    int bmpWidth = channel_s->width;
//    int bmpHeight = channel_s->height;
//    int bmpStep = channel_s->widthStep;
//
//    for(j=0;j<bmpHeight;j++)//y
//    {
//        for(i=0;i<bmpWidth;i++)//x
//        {
//            channel_s->imageData[j*bmpStep+i] -
//        }
//    }



//
//    CvMoments moments;
//    CvHuMoments hu;
//    cvMoments(channel_h,&moments,0);
//    cvGetHuMoments(&moments, &hu);





//    Moments m = moments((uchar*)channel_h->imageData);
//    cvMoments()
//    0.00432353 4.09781e-10 3.51024e-13 7.65516e-12 -9.62173e-25 -4.51895e-17 -1.25118e-23
//    0.00144638 5.4656e-11 4.14596e-14 2.39258e-13 -1.1285e-26 7.43172e-19 -2.09877e-26
//    0.00144187 6.58407e-11 5.8827e-14 3.43219e-13 -2.74641e-26 4.2163e-19 -4.03007e-26
    //ustom hu
//    HuMoment(&gray);
//    Moments m = moments(gray);
    //opencv hu

//    for(int i = 0; i< gray.rows;i++){
//        for(int j = 0;j<gray.cols;j++){
//            cout << int(gray.at<uchar>(i,j)) << endl;
//        }
//    }
    return 0;
}

void calcCiCS(vector<vector<double *>> * datas,IplImage *channel_h,IplImage *channel_s,IplImage *channel_v){
    vector<double> Hus_h = HuMoment(channel_h);
    vector<double> Hus_s = HuMoment(channel_s);
    vector<double> Hus_v = HuMoment(channel_v);
    int i,j;

    int bmpWidth = channel_h->width;
    int bmpHeight = channel_h->height;
    int bmpStep = channel_h->widthStep;
//    cout << "complete" << Hus_h[0] <<endl;

    double sum_h = Hus_h[1] + Hus_h[2] + Hus_h[3] + Hus_h[4];
    double sum_s = Hus_s[1] + Hus_s[2] + Hus_s[3] + Hus_s[4];
    double sum_v = Hus_v[1] + Hus_v[2] + Hus_v[3] + Hus_v[4];

    for(j=0;j<bmpHeight;j++)//y
    {
        vector<double *> data;
        for(i=0;i<bmpWidth;i++)//x
        {
            cout <<  int(channel_h->imageData[j*bmpStep+i]) << " ";
            channel_h->imageData[j*bmpStep+i] = channel_h->imageData[j*bmpStep+i] - sum_h;
            channel_s->imageData[j*bmpStep+i] = channel_s->imageData[j*bmpStep+i] - sum_s;
            channel_v->imageData[j*bmpStep+i] = channel_v->imageData[j*bmpStep+i] - sum_v;
            cout <<  int(channel_h->imageData[j*bmpStep+i]);
            //hen,calc cics1 cics2
            //cics1 = cicsH/cicsS
            //cics2 = cicsS/cicsV
            double d[2];
            d[0] = channel_h->imageData[j*bmpStep+i]/channel_s->imageData[j*bmpStep+i];
            d[1] = channel_s->imageData[j*bmpStep+i]/channel_v->imageData[j*bmpStep+i];
            cout << " " << d[0] << " " << d[1] << endl;
            data.push_back(d);
        }

        datas->push_back(data);
    }
//    cout << "complete" <<endl;
}

//#################################################################################//
double * getHu(Moments m){
    double M[7] = {0};        //HU不变矩
    //y20 - y02
    double t1=(m.nu20-m.nu02);
    // y30 - 3 * y12
    double t2=(m.nu30-3*m.nu12);
    // 3 * y21 - y03
    double t3=(3*m.nu21-m.nu03);
    // y30 + y12
    double t4=(m.nu30+m.nu12);
    // y21 + y03
    double t5=(m.nu21+m.nu03);
    M[0]=m.nu20+m.nu02;
    M[1]=t1*t1+4*m.nu11*m.nu11;
    M[2]=t2*t2+t3*t3;
    M[3]=t4*t4+t5*t5;
    M[4]=t2*t4*(t4*t4-3*t5*t5)+t3*t5*(3*t4*t4-t5*t5);
    M[5]=t1*(t4*t4-t5*t5)+4*m.nu11*t4*t5;
    M[6]=t3*t4*(t4*t4-3*t5*t5)-t2*t5*(3*t4*t4-t5*t5);
    return M;
}
bool HuMoment(Mat* img)
{
    double M[7] = {0};
    int bmpWidth = img->rows;
    int bmpHeight = img->cols;
//    int bmpStep = img->widthStep;
    int bmpChannels = img->channels();
//    uchar* pBmpBuf = (uchar*)img->imageData;

    double m00=0,m11=0,m20=0,m02=0,m30=0,m03=0,m12=0,m21=0;  //中心矩
    double x0=0,y0=0;    //计算中心距时所使用的临时变量（x-x'）
    double u20=0,u02=0,u11=0,u30=0,u03=0,u12=0,u21=0;//规范化后的中心矩
    //double M[7];    //HU不变矩
    double t1=0,t2=0,t3=0,t4=0,t5=0;//临时变量，
    //double Center_x=0,Center_y=0;//重心
    int Center_x=0,Center_y=0;//重心
    int i,j;            //循环变量

    //  获得图像的区域重心(普通矩)
    double s10=0,s01=0,s00=0;  //0阶矩和1阶矩
    for(j=0;j<bmpHeight;j++)//y
    {
        for(i=0;i<bmpWidth;i++)//x
        {
//            ;
//            s10+=i*pBmpBuf[j*bmpStep+i];
//            s01+=j*pBmpBuf[j*bmpStep+i];
//            s00+=pBmpBuf[j*bmpStep+i];
            s10+=i*int(img->at<uchar>(i,j));
            s01+=j*int(img->at<uchar>(i,j));
            s00+=int(img->at<uchar>(i,j));
        }
    }
    Center_x=(int)(s10/s00+0.5);
    Center_y=(int)(s01/s00+0.5);

    //  计算二阶、三阶矩(中心矩)
    m00=s00;
    for(j=0;j<bmpHeight;j++)
    {
        for(i=0;i<bmpWidth;i++)//x
        {
            x0=(i-Center_x);
            y0=(j-Center_y);
            m11+=x0*y0*int(img->at<uchar>(i,j));
            m20+=x0*x0*int(img->at<uchar>(i,j));
            m02+=y0*y0*int(img->at<uchar>(i,j));
            m03+=y0*y0*y0*int(img->at<uchar>(i,j));
            m30+=x0*x0*x0*int(img->at<uchar>(i,j));
            m12+=x0*y0*y0*int(img->at<uchar>(i,j));
            m21+=x0*x0*y0*int(img->at<uchar>(i,j));
        }
    }

    //  计算规范化后的中心矩: mij/pow(m00,((i+j+2)/2)
    u20=m20/pow(m00,2);
    u02=m02/pow(m00,2);
    u11=m11/pow(m00,2);
    u30=m30/pow(m00,2.5);
    u03=m03/pow(m00,2.5);
    u12=m12/pow(m00,2.5);
    u21=m21/pow(m00,2.5);

    //  计算中间变量
    //y20 - y02
    t1=(u20-u02);
    // y30 - 3 * y12
    t2=(u30-3*u12);
    // 3 * y21 - y03
    t3=(3*u21-u03);
    // y30 + y12
    t4=(u30+u12);
    // y21 + y03
    t5=(u21+u03);

    //  计算不变矩
    M[0]=u20+u02;
    M[1]=t1*t1+4*u11*u11;
    M[2]=t2*t2+t3*t3;
    M[3]=t4*t4+t5*t5;
    M[4]=t2*t4*(t4*t4-3*t5*t5)+t3*t5*(3*t4*t4-t5*t5);
    M[5]=t1*(t4*t4-t5*t5)+4*u11*t4*t5;
    M[6]=t3*t4*(t4*t4-3*t5*t5)-t2*t5*(3*t4*t4-t5*t5);


    cout << M[0] << " " << M[1]<< " " << M[2]<< " " << M[3]
            << " " << M[4]<< " " << M[5]<< " " << M[6]<< endl;
    return true;
}
vector<double> HuMoment(IplImage* img)
{
    double M[7] = {0};
    int bmpWidth = img->width;
    int bmpHeight = img->height;
    int bmpStep = img->widthStep;
    int bmpChannels = img->nChannels;
    uchar*pBmpBuf = (uchar*)img->imageData;

    double m00=0,m11=0,m20=0,m02=0,m30=0,m03=0,m12=0,m21=0;  //中心矩
    double x0=0,y0=0;    //计算中心距时所使用的临时变量（x-x'）
    double u20=0,u02=0,u11=0,u30=0,u03=0,u12=0,u21=0;//规范化后的中心矩
    //double M[7];    //HU不变矩
    double t1=0,t2=0,t3=0,t4=0,t5=0;//临时变量，
    //double Center_x=0,Center_y=0;//重心
    int Center_x=0,Center_y=0;//重心
    int i,j;            //循环变量

    //  获得图像的区域重心(普通矩)
    double s10=0,s01=0,s00=0;  //0阶矩和1阶矩
    for(j=0;j<bmpHeight;j++)//y
    {
        for(i=0;i<bmpWidth;i++)//x
        {
            s10+=i*pBmpBuf[j*bmpStep+i];
            s01+=j*pBmpBuf[j*bmpStep+i];
            s00+=pBmpBuf[j*bmpStep+i];
        }
    }
    Center_x=(int)(s10/s00+0.5);
    Center_y=(int)(s01/s00+0.5);

    //  计算二阶、三阶矩(中心矩)
    m00=s00;
    for(j=0;j<bmpHeight;j++)
    {
        for(i=0;i<bmpWidth;i++)//x
        {
            x0=(i-Center_x);
            y0=(j-Center_y);
            m11+=x0*y0*pBmpBuf[j*bmpStep+i];
            m20+=x0*x0*pBmpBuf[j*bmpStep+i];
            m02+=y0*y0*pBmpBuf[j*bmpStep+i];
            m03+=y0*y0*y0*pBmpBuf[j*bmpStep+i];
            m30+=x0*x0*x0*pBmpBuf[j*bmpStep+i];
            m12+=x0*y0*y0*pBmpBuf[j*bmpStep+i];
            m21+=x0*x0*y0*pBmpBuf[j*bmpStep+i];
        }
    }

    //  计算规范化后的中心矩: mij/pow(m00,((i+j+2)/2)
    u20=m20/pow(m00,2);
    u02=m02/pow(m00,2);
    u11=m11/pow(m00,2);
    u30=m30/pow(m00,2.5);
    u03=m03/pow(m00,2.5);
    u12=m12/pow(m00,2.5);
    u21=m21/pow(m00,2.5);

    //  计算中间变量
    t1=(u20-u02);
    t2=(u30-3*u12);
    t3=(3*u21-u03);
    t4=(u30+u12);
    t5=(u21+u03);

    //  计算不变矩
    vector<double> ret;
    M[0]=u20+u02;
    M[1]=t1*t1+4*u11*u11;
    M[2]=t2*t2+t3*t3;
    M[3]=t4*t4+t5*t5;
    M[4]=t2*t4*(t4*t4-3*t5*t5)+t3*t5*(3*t4*t4-t5*t5);
    M[5]=t1*(t4*t4-t5*t5)+4*u11*t4*t5;
    M[6]=t3*t4*(t4*t4-3*t5*t5)-t2*t5*(3*t4*t4-t5*t5);
    cout << M[0] << " " << M[1]<< " " << M[2]<< " " << M[3]
         << " " << M[4]<< " " << M[5]<< " " << M[6]<< endl;

    ret.push_back(M[0]);
    ret.push_back(M[1]);
    ret.push_back(M[2]);
    ret.push_back(M[3]);
    ret.push_back(M[4]);
    ret.push_back(M[5]);
    ret.push_back(M[6]);
    return ret;
}
