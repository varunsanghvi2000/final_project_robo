
#include <PID_v1.h>


#define analogInPin1 A0  // Analog input pin that the potentiometer is attached to
#define analogInPin2 A1  // Analog input pin that the potentiometer is attached to
#define analogInPin3 A2  // Analog input pin that the potentiometer is attached to

#define analogInPinop1 A3  // Analog input pin that the potentiometer is attached to
#define analogInPinop2 A4  // Analog input pin that the potentiometer is attached to
#define analogInPinop3 A5  // Analog input pin that the potentiometer is attached to

#define arm1_left 0
#define arm1_right 1
#define arm2_left 2
#define arm2_right 3
#define arm3_left 4
#define arm3_right 5

#define arm1_pwm 9 
#define arm2_pwm 10
#define arm3_pwm 14

double kp = 2 , ki = 1 , kd = 0.011;             // modify for optimal performance
double input1 = 0, output1 = 0, setpoint1 = 0, offset1=-343;
double input2 = 0, output2 = 0, setpoint2 = 0, offset2=68;
double input3 = 0, output3 = 0, setpoint3 = 0, offset3=-75;

PID PID1(&input1, &output1, &setpoint1, kp, ki, kd, DIRECT);  // if motor will only run at full speed try 'REVERSE' instead of 'DIRECT'
PID PID2(&input2, &output2, &setpoint2, kp, ki, kd, DIRECT);  // if motor will only run at full speed try 'REVERSE' instead of 'DIRECT'
PID PID3(&input3, &output3, &setpoint3, kp, ki, kd, DIRECT);  // if motor will only run at full speed try 'REVERSE' instead of 'DIRECT'

//const int analogOutPin = 9; // Analog output pin that the LED is attached to

int sl_arm1 = 0;        // value read from the pot
int sl_arm2 = 0;        // value read from the pot
int sl_arm3 = 0;        // value read from the pot
int man_arm1= 0;        // value read from the pot
int man_arm2 = 0;        // value read from the pot
int man_arm3 = 0;        // value read from the pot


//int arm1_pwm = 0;        // value output to the PWM (analog out)
//int arm2_pwm = 0;        // value output to the PWM (analog out)
//int arm3_pwm = 0;        // value output to the PWM (analog out)


// robot geometry
 // (look at pics above for explanation)
 const float e = 115.0;     // end effector
 const float f = 457.3;     // base
 const float re = 232.0;
 const float rf = 112.0;
 
 // trigonometric constants
 const float sqrt3 = sqrt(3.0);
 const float pi = 3.141592653;    // PI
 const float sin120 = sqrt3/2.0;   
 const float cos120 = -0.5;        
 const float tan60 = sqrt3;
 const float sin30 = 0.5;
 const float tan30 = 1/sqrt3;
 
 // forward kinematics: (theta1, theta2, theta3) -> (x0, y0, z0)
 // returned status: 0=OK, -1=non-existing position
 int delta_calcForward(float theta1, float theta2, float theta3, float &x0, float &y0, float &z0) {
     float t = (f-e)*tan30/2;
     float dtr = pi/(float)180.0;
 
     theta1 *= dtr;
     theta2 *= dtr;
     theta3 *= dtr;
 
     float y1 = -(t + rf*cos(theta1));
     float z1 = -rf*sin(theta1);
 
     float y2 = (t + rf*cos(theta2))*sin30;
     float x2 = y2*tan60;
     float z2 = -rf*sin(theta2);
 
     float y3 = (t + rf*cos(theta3))*sin30;
     float x3 = -y3*tan60;
     float z3 = -rf*sin(theta3);
 
     float dnm = (y2-y1)*x3-(y3-y1)*x2;
 
     float w1 = y1*y1 + z1*z1;
     float w2 = x2*x2 + y2*y2 + z2*z2;
     float w3 = x3*x3 + y3*y3 + z3*z3;
     
     // x = (a1*z + b1)/dnm
     float a1 = (z2-z1)*(y3-y1)-(z3-z1)*(y2-y1);
     float b1 = -((w2-w1)*(y3-y1)-(w3-w1)*(y2-y1))/2.0;
 
     // y = (a2*z + b2)/dnm;
     float a2 = -(z2-z1)*x3+(z3-z1)*x2;
     float b2 = ((w2-w1)*x3 - (w3-w1)*x2)/2.0;
 
     // a*z^2 + b*z + c = 0
     float a = a1*a1 + a2*a2 + dnm*dnm;
     float b = 2*(a1*b1 + a2*(b2-y1*dnm) - z1*dnm*dnm);
     float c = (b2-y1*dnm)*(b2-y1*dnm) + b1*b1 + dnm*dnm*(z1*z1 - re*re);
  
     // discriminant
     float d = b*b - (float)4.0*a*c;
     if (d < 0) return -1; // non-existing point
 
     z0 = -(float)0.5*(b+sqrt(d))/a;
     x0 = (a1*z0 + b1)/dnm;
     y0 = (a2*z0 + b2)/dnm;
     return 0;
 }
 
 // inverse kinematics
 // helper functions, calculates angle theta1 (for YZ-pane)
 int delta_calcAngleYZ(float x0, float y0, float z0, float &theta) {
     float y1 = -0.5 * 0.57735 * f; // f/2 * tg 30
     y0 -= 0.5 * 0.57735    * e;    // shift center to edge
     // z = a + b*y
     float a = (x0*x0 + y0*y0 + z0*z0 +rf*rf - re*re - y1*y1)/(2*z0);
     float b = (y1-y0)/z0;
     // discriminant
     float d = -(a+b*y1)*(a+b*y1)+rf*(b*b*rf+rf); 
     if (d < 0) return -1; // non-existing point
     float yj = (y1 - a*b - sqrt(d))/(b*b + 1); // choosing outer point
     float zj = a + b*yj;
     theta = 180.0*atan(-zj/(y1 - yj))/pi + ((yj>y1)?180.0:0.0);
     return 0;
 }
 
 // inverse kinematics: (x0, y0, z0) -> (theta1, theta2, theta3)
 // returned status: 0=OK, -1=non-existing position
 int delta_calcInverse(float x0, float y0, float z0, float &theta1, float &theta2, float &theta3) {
     theta1 = theta2 = theta3 = 0;
     int status = delta_calcAngleYZ(x0, y0, z0, theta1);
     if (status == 0) status = delta_calcAngleYZ(x0*cos120 + y0*sin120, y0*cos120-x0*sin120, z0, theta2);  // rotate coords to +120 deg
     if (status == 0) status = delta_calcAngleYZ(x0*cos120 - y0*sin120, y0*cos120+x0*sin120, z0, theta3);  // rotate coords to -120 deg
     return status;
 }


void pid_setup(){
  PID1.SetMode(AUTOMATIC);
  PID1.SetSampleTime(1);
  PID1.SetOutputLimits(-255, 255);
  PID2.SetMode(AUTOMATIC);
  PID2.SetSampleTime(1);
  PID2.SetOutputLimits(-255, 255);
  PID3.SetMode(AUTOMATIC);
  PID3.SetSampleTime(1);
  PID3.SetOutputLimits(-255, 255);
}

void drivemotor(int out1,int out2,int out3)
{
  if (out1 > 0) {
    digitalWrite(arm2_left,HIGH);
    digitalWrite(arm2_right,LOW);
    analogWrite(arm2_pwm, out1);
  }
  else {
    digitalWrite(arm2_left,LOW);
    digitalWrite(arm2_right,HIGH);
    analogWrite(arm2_pwm, abs(out1));
  }
  if (out2 > 0) {
    digitalWrite(arm1_left,LOW);
    digitalWrite(arm1_right,HIGH);
    analogWrite(arm1_pwm, out2);
  }
  else {
    digitalWrite(arm1_left,HIGH);
    digitalWrite(arm1_right,LOW);
    analogWrite(arm1_pwm, abs(out2));
  }
  if (out3 > 0) {
    digitalWrite(arm3_left,LOW);
    digitalWrite(arm3_right,HIGH);
    analogWrite(arm3_pwm, out3);
  }
  else {
    digitalWrite(arm3_left,HIGH);
    digitalWrite(arm3_right,LOW);
    analogWrite(arm3_pwm, abs(out3));
  }
  
}


void setup() {
  // initialize serial communications at 9600 bps:
  Serial.begin(9600);
  pid_setup();
  pinSetup();
}

void loop() {
  setpoint1 = analogRead(analogInPin1);                       // modify to fit motor and encoder characteristics, potmeter connected to A0
  delay(5);
  setpoint2 = analogRead(analogInPin2);                       // modify to fit motor and encoder characteristics, potmeter connected to A0
  delay(5);
  setpoint3 = analogRead(analogInPin3);                       // modify to fit motor and encoder characteristics, potmeter connected to A0
  delay(5);

//  setpoint1 = 816;                       // modify to fit motor and encoder characteristics, potmeter connected to A0
//  //delay(5);
//  setpoint2 = 488;                       // modify to fit motor and encoder characteristics, potmeter connected to A0
//  //delay(5);
//  setpoint3 = 557;                       // modify to fit motor and encoder characteristics, potmeter connected to A0
//  //delay(5);


  input1 = analogRead(analogInPinop1) - offset1;              // data from encoder
  delay(5);
  input2 = analogRead(analogInPinop2) - offset2;
  delay(5);
  input3 = analogRead(analogInPinop3) - offset3;
  delay(5);
  Serial.print("sensor 1 ");
  Serial.print(setpoint1);
  Serial.print("sensor op 1 ");
  Serial.print(input1);
  Serial.print("sensor 2 ");
  Serial.print(setpoint2);
  Serial.print("sensor op 2 ");
  Serial.print(input2);
  Serial.print("sensor 3 ");
  Serial.print(setpoint3);
  Serial.print("sensor op 3 ");
  Serial.println(input3);
//    
  PID1.Compute();                                             // calculate new output
  PID2.Compute();                                             // calculate new output
  PID3.Compute();                                             // calculate new output
  drivemotor(output1,output2,output3);
//  Serial.print("output1");
//  Serial.println(output1);
//    digitalWrite(arm1_left,LOW);
//    digitalWrite(arm1_right,HIGH);
//    analogWrite(arm1_pwm,255);
//    digitalWrite(arm2_left,HIGH);
//    digitalWrite(arm2_right,LOW);
//    analogWrite(arm2_pwm,255);
//    digitalWrite(arm3_left,LOW);
//    digitalWrite(arm3_right,HIGH);
//    analogWrite(arm3_pwm,255);
delay(50); 
  
}


void pinSetup()
{
pinMode(arm1_left,OUTPUT);
pinMode(arm1_right,OUTPUT);
pinMode(arm2_left,OUTPUT);
pinMode(arm2_right,OUTPUT);
pinMode(arm3_left,OUTPUT);
pinMode(arm3_right,OUTPUT);
//pinMode(arm1_pwm,OUTPUT);
//pinMode(arm2_pwm,OUTPUT);
//pinMode(arm3_pwm,OUTPUT);
digitalWrite(arm1_left,LOW);
digitalWrite(arm1_right,LOW);
digitalWrite(arm2_left,LOW);
digitalWrite(arm2_right,LOW);
digitalWrite(arm3_left,LOW);
digitalWrite(arm3_right,LOW);
//digitalWrite(arm1_pwm,LOW);
//digitalWrite(arm2_pwm,LOW);
//digitalWrite(arm3_pwm,LOW);


}


