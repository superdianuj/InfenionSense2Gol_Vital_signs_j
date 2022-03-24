void setup() {
  // put your setup code here, to run once:

  Serial.begin(9600);
  pinMode(14,INPUT);
  pinMode(12,INPUT);

}

void loop() {
  // put your main code here, to run repeatedly:

  if((digitalRead(10)==1)||(digitalRead(11)==1))
  {
    Serial.println('!');
    }
    else
    {
      Serial.println(analogRead(A0));
      }


      delay(1);

}
