
// Read read data poits from a file "FileName"

int ReadData(char* FileName, float x[], float xscale, float y[], float yscale )
{
 
 ifstream in;
 in.open(FileName);
 int nlines = 0; 

while (1) {
      in >> x[nlines] >> y[nlines];
      if (!in.good()) break;
      x[nlines]=xscale*x[nlines];
      y[nlines]=yscale*y[nlines];
      nlines++;
   }
 cout<<"*** In total "<<nlines<<" data points taken from "<<FileName<<" and rescaled"<<endl;
 cout<<"    First point: x="<<x[0]<<"  y="<<y[0]<<endl;
 cout<<"    Last point:  x="<<x[nlines-1]<<"  y="<<y[nlines-1]<<endl;
 return nlines;
}
