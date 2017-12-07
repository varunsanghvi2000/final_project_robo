// include libraries
#include <WiFiClient.h>
#include <ESP8266WebServer.h>
#include <SPI.h>
#include <ESP8266WiFi.h>
#include <WiFiUdp.h>

// wifi variables
const char* ssid = "NETGEAR68";
const char* password = "JOLLYMINT969";

unsigned int remoteport=1390 ;      // local port to listen on
unsigned int localPort =  1808 ;      // local port to listen on

IPAddress remoteip(192, 168, 1, 64);
WiFiUDP Udp;


int value,r,g,b;
int leftint=0;
int rightint=0;
int leftdata=0;
int rightdata=0;
char packetBuffer[255]; //buffer to hold incoming packet
char command_init[1];
char  ReplyBuffer[] = "141414141";       // a string to send back
char *pt;
char *sec;
float acc=0;
float til=0;
float ltilt=0;
float rtilt=0;
int motleft=0;
int motright=0;
int comman=0;

void setup() {
  setupSerial();
  pinsetup();
  setupudp();
  motorstop();
}

void loop() {
  // put your main code here, to run repeatedly:
readdata(leftdata,rightdata);
}



void readdata(int leftdata,int rightdata)
{
    int packetSize = Udp.parsePacket();
    if (packetSize) {
    int len = Udp.read(packetBuffer, 255);
    if (len > 0) {
      packetBuffer[len] = 0;
    }
    //command_init[0]=packetBuffer[0];
    pt = strtok (packetBuffer,":");
    sec=pt;
    Serial.println("command ");
    Serial.println(sec);
    comman=atoi(sec);
    Serial.println("com ");
    Serial.println(comman);
    pt = strtok (NULL,":");
    Serial.println("value ");
    //Serial.println(pt);
    value=atoi(pt);
    Serial.println(value);

    if(comman==1)
    {
      Serial.println("comman ");
      Serial.println(comman);
      Serial.println("here1 ");
      Serial.println(value);
      if(value==0)
        motorstop();
      if(value==1)
        motorforward();
      if(value==2)
        motorreverse();
    }else if(comman==2)
    {
            Serial.println("comman ");
      Serial.println(comman);

      
       acc=(float)value/10;
       Serial.println("acc ");
      Serial.println(acc);
    }else if(comman==3)
    {
      Serial.println("comman ");
      Serial.println(comman);
      Serial.println("value");
      Serial.println(value);
       til=(float)value/100;
       Serial.println("value/100");
       Serial.println(til);
       if(til<0.9)
       {
       ltilt=1;
       rtilt=til+0.1; 
       Serial.println("r ltilt ");
      Serial.println(ltilt);
      Serial.println("r rtilt");
      Serial.println(rtilt);
       }else if(til>0.9)
       {
       ltilt=1-abs(til-0.9)+0.1;
       rtilt=1; 
       Serial.println("l ltilt ");
       Serial.println(ltilt);
       Serial.println("l rtilt");
       Serial.println(rtilt);
       }else{
       ltilt=1;
       rtilt=1; 
       //       Serial.println("ltilt ");
     // Serial.println(ltilt);
   //   Serial.println("rtilt");
  //    Serial.println(rtilt);
       }
    }
    
  }
  
void pinsetup()
{
  pinMode(motor1Enable, OUTPUT);
  pinMode(motor1Forward, OUTPUT);
  pinMode(motor1Back, OUTPUT);
  pinMode(motor2Enable, OUTPUT);
  pinMode(motor2Forward, OUTPUT);
  pinMode(motor2Back, OUTPUT);
  pinMode(trigPin, INPUT);
}

void setupudp()
{
  if(digitalRead(trigPin)==HIGH)
  {
  WiFi.begin(ssid0, password0);// Connect to WiFi 
  WiFi.config(IPAddress(192, 168, 1, 64),
            IPAddress(192, 168, 1, 1),
              IPAddress(255, 255, 255, 0));
  }else
  {
  WiFi.begin(ssid0, password0);// Connect to WiFi 
  WiFi.config(IPAddress(192, 168, 1, 64),
  IPAddress(192, 168, 1, 1),
  IPAddress(255, 255, 255, 0));
  }
  // while wifi not connected yet, print '.'
  // then after it connected, get out of the loop
  while (WiFi.status() != WL_CONNECTED) {
     delay(500);
     Serial.print(".");
  }
  //print a new line, then print WiFi connected and the IP address
  Serial.println("");
  Serial.println("WiFi connected");
  // Print the IP address
  Serial.println(WiFi.localIP());
  Udp.begin(localPort);

}

void setupSerial()
{
    Serial.begin(9600);
  
  while (!Serial) {
    ; // wait for serial port to connect. Needed for native USB port only
  }
}

void motorstop()
{
    //motors stop
    //send_drive_command
//    digitalWrite(motor1Forward, LOW);
//  digitalWrite(motor1Back, LOW);
//  digitalWrite(motor2Forward, LOW);
//  digitalWrite(motor2Back, LOW);
}


void motorforward()
{
    //motors forward
//  digitalWrite(motor1Forward, HIGH);
//  digitalWrite(motor1Back, LOW);
//  digitalWrite(motor2Forward, HIGH);
//  digitalWrite(motor2Back, LOW);
}


void motorreverse()
{
    //motors forward
//  digitalWrite(motor1Forward, LOW);
//  digitalWrite(motor1Back, HIGH);
//  digitalWrite(motor2Forward, LOW);
//  digitalWrite(motor2Back, HIGH);
}
