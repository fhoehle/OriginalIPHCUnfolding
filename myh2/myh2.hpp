#ifndef MYH2_HPP
#define MYH2_HPP

#include <iostream>
#include <cmath>
#include <algorithm>

#include <TH2F.h>
#include <TString.h>

using namespace std;

// Yes, all of these methods use proper camelCase instead of ROOT's stupid capitalization.
// Why? Well, besides looking better it should help avoid confusion when porting code from
// TH2 to MyH2, because compiler errors will be triggered when you weren't careful enough.

/// 2 dimensional histogram which allows to have different y-bin-edges for different x-bins;
/// so for example a histogram can have y bin edges 0,3,5 in the first x bin but
///                                                 0,4,5 in the second x bin and so on.
/// There are no overflow bins! Events which are too far out will be filled into the outermost bins.
template <int nXBin, int nYBin>
class MyH2
{
  public:
    
    /// constructor
    MyH2(const double* xBinEdges_, const double yBinEdges_[nXBin][nYBin+1])
      : name("")
    {
      init(xBinEdges_, yBinEdges_);
    }
    
    /// constructor with name
    MyH2(const char* name_, const double* xBinEdges_, const double yBinEdges_[nXBin][nYBin+1])
      : name(name_)
    {
      init(xBinEdges_, yBinEdges_);
    }
    
    /// constructor with two names, for compat with stupid TH2 - we'll ignore the first name
    MyH2(const char* , const char* name_, const double* xBinEdges_, const double yBinEdges_[nXBin][nYBin+1])
      : name(name_)
    {
      init(xBinEdges_, yBinEdges_);
    }
    
    
    /// read values and errors from a TH2 - this is how you convert a TH2 into a MyH2
    void readTH2(TH2* h)
    {
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          values[x][y] = h->GetBinContent(x+1,y+1);
          setBinError(x, y, h->GetBinError(x+1,y+1));
        }
      }      
    }
    
    /// generate TH2F from this, for example for saving to a file or doing anything that
    ///  doesn't depend on bin edges or errors
    TH2F* toTH2F(const char* name_)
    {
      TH2F* h = new TH2F(name_, name_, nXBin, 0.5, nXBin+0.5, nYBin, 0.5, nYBin+0.5);
      
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          h->SetBinContent(x+1, y+1, values[x][y]);
          h->SetBinError(x+1, y+1, getBinError(x,y));
        }
      }
      return h;
    }
    
    /// generate TH2F from this, for example for saving to a file or doing anything that
    ///  doesn't depend on bin edges or errors
    TH2F* toTH2F()
    {
      return toTH2F(name);
    }
    
    /// generate histogram that looks approx. like this MyH2 would if it was plottable
    TH2F* toVisTH2F(const char* name_)
    {
      const double miny = yBinEdges[0][0];
      const double maxy = yBinEdges[0][nYBin];
      
      TH2F* h = new TH2F(name_, name_, nXBin, xBinEdges, 1000, miny, maxy);
      
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<1000; y++)
        {
          const int xbin = x;
          const int ybin = findYBin(xbin, h->GetYaxis()->GetBinCenter(y+1));
          
          
          h->SetBinContent(x+1, y+1, values[xbin][ybin]);
          h->SetBinError(x+1, y+1, getBinError(xbin,ybin));
        }
      }
      return h;
    }
    
    /// generate histogram that looks approx. like this MyH2 would if it was plottable
    TH2F* toVisTH2F()
    {
      return toVisTH2F(name);
    }
    
    /// add entry, overflow events are taken into outermost bins
    void fill(const double x, const double y, const double weight)
    {
      // find x index
      const int xbin = findXBin(x);
      
      // find y index
      const int ybin = findYBin(xbin, y);
      
      // increment
      values[xbin][ybin] += weight;
      w2matrix[xbin][ybin] += weight*weight;
 
    }
    
    /*// old code that didn't do a binary search, should be slower than what we have now
    /// get bin x index for specific x value
    int findXBin(double x)
    {
      int xbin = 0;
      for(int i=0; i<nXBin; i++)
      {
        if(x>xBinEdges[i+1])
          continue;
        xbin = i;
        break;
      }
      if(x>xBinEdges[nXBin])
        xbin = nXBin-1;
      
      return xbin;
    }
    
    /// get bin y index for y value and bin x index
    int findYBin(int xbin, double y)
    {      
      int ybin = 0;
      for(int i=0; i<nYBin; i++)
      {
        if(y>yBinEdges[xbin][i+1])
          continue;
        ybin = i;
        break;
      }
      if(y>yBinEdges[xbin][nYBin])
        ybin = nYBin-1;
      
      return ybin;
    }*/
    
    /// get bin x index for specific x value
    int findXBin(const double x)
    {
      const int xbin = int(lower_bound(xBinEdges, xBinEdges + nXBin+1, x) - xBinEdges -1);

      if(xbin==nXBin)
        return nXBin-1;
      else if(xbin==-1) 
        return 0;
      else
        return xbin;
    }
    
    /// get bin y index for y value and bin x index
    int findYBin(const int xbin, const double y)
    {      
      const int ybin = int(lower_bound(yBinEdges[xbin], yBinEdges[xbin] + nYBin+1, y) - yBinEdges[xbin] -1);
      
      if(ybin==nYBin)
        return nYBin-1;
      else if(ybin==-1) 
        return 0;
      else
        return ybin;
    }
    
    int findYBin(double xbin, double y); // hack: undefined functions so nobody can accidentally
    int findYBin(float xbin, double y); // pass a real x value instead of a bin number. won't link then.
    
    
    /// bin content getter
    double getBinContent(int x, int y)
    {
      return values[x][y];
    }
    
    /// bin error getter
    double getBinError(int x, int y)
    {
      return sqrt(w2matrix[x][y]);
    }
    
    /// name getter
    TString getName()
    {
      return name;
    }
    
    /// bin content setter
    void setBinContent(int x, int y, double content)
    {
      values[x][y] = content;
    }
    
    /// bin error setter
    void setBinError(int x, int y, double error)
    {
      w2matrix[x][y] = error*error;
    }
    
    /// name setter
    void setName(const char* name_)
    {
      name = name_;
    }
    
    
    /// gets center of bin on the x axis, zero-based index
    double getXBinCenter(int x)
    {
      return (xBinEdges[x+1]+xBinEdges[x])/2;
    }
    
    /// integral
    double integral(int binx1, int binx2, int biny1, int biny2)
    {
      double result = 0;
      for(int x=binx1; x<=binx2; x++)
      {
        for(int y=biny1; y<=biny2; y++)
        {
          result += values[x][y];
        }
      }
      return result;
    }
    
    /// total integral
    double integral()
    {
      return integral(0, nXBin-1, 0, nYBin-1);
    }
    
    /// add other MyH2, with scale factor for the other histo
    void add(MyH2* otherh, double scale_=1)
    {
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          values[x][y] += otherh->getBinContent(x,y) * scale_;
          const double error = otherh->getBinError(x,y);
          w2matrix[x][y] += error * error * scale_ * scale_;
        }
      }
      
    }
    
    /// add TH2, with scale factor for the other histo
    void add(TH2* otherh, double scale_=1)
    {
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          values[x][y] += otherh->GetBinContent(x+1,y+1) * scale_;
          const double error = otherh->GetBinError(x+1,y+1);
          w2matrix[x][y] += error * error * scale_ * scale_;
        }
      }
      
    }
    
    /// divide by other MyH2
    void divide(MyH2* otherh)
    {
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          const double content = values[x][y];
          const double otherContent = otherh->getBinContent(x,y);
          const double error2 = w2matrix[x][y];
          const double otherError2 = otherh->getBinError(x,y) * otherh->getBinError(x,y);
          
          if(otherContent)
            values[x][y] /= otherContent;  
          else
            values[x][y] = 0;
          
          if(!otherContent)
            w2matrix[x][y] = 0;
          else
          {            
            w2matrix[x][y]  = (error2 * otherContent * otherContent);
            w2matrix[x][y] += (otherError2 * content * content);
            w2matrix[x][y] /= (otherContent * otherContent * otherContent * otherContent);
          }
        }
      }
      
    }
    
    /// multiply a global factor to histogram
    void scale(double scale_)
    {
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          values[x][y] *= scale_;
          
          w2matrix[x][y] *= scale_ * scale_;
        }
      }
    }
    
    /// make a new clone of this MyH2
    MyH2<nXBin, nYBin>* clone()
    {
      MyH2<nXBin, nYBin>* h = new MyH2<nXBin, nYBin>(name, xBinEdges, yBinEdges);
      
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          h->setBinContent(x,y, values[x][y]);
          h->setBinError(x,y, getBinError(x,y));
        }
      }
      
      h->setName(name);
      
      return h;
    }
    
    /// emulate TH2's Write by converting to TH2, calling Write, deleting the new object again
    void write()
    {
      TH2F* h = toTH2F();
      h->Write();
      delete h;
    }
    
    /// debug output
    void printXEdges()
    {
      for(int x=0; x<=nXBin; x++)
      {
        cout << xBinEdges[x] << " | ";
      }
      cout << endl;
    }
  
  
    /// debug output
    void printYEdges()
    {
      for(int x=0; x<nXBin; x++)
      {
        cout << "xBin " << x << "\n\t";
        for(int y=0; y<=nYBin; y++)
        {
          cout << yBinEdges[x][y] << " | ";
        }
        cout << endl;
      }
    }
    
    /// debug output
    void printValues()
    {
      for(int x=0; x<nXBin; x++)
      {
        cout << "xBin " << x << "\n\t";
        for(int y=0; y<nYBin; y++)
        {
          cout << values[x][y] << " | ";
        }
        cout << endl;
      }
    }
  
  private:
    
    // initialize array edges - needs to be a separate function because constructors can't call each other
    // in old versions of c++
    void init(const double* xBinEdges_, const double yBinEdges_[nXBin][nYBin+1])
    {
      // get x edges
      for(int i=0; i<=nXBin; i++)
      {
        xBinEdges[i] = xBinEdges_[i];
      }
      
      // get y edges in all variants
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<=nYBin; y++)
        {
          yBinEdges[x][y] = yBinEdges_[x][y];
        }
      }
      
      // init to 0
      for(int x=0; x<nXBin; x++)
      {
        for(int y=0; y<nYBin; y++)
        {
          values[x][y] = 0;
          w2matrix[x][y] = 0;
        }
      }      
    }
    
    
    
    double xBinEdges[nXBin+1];
    double yBinEdges[nXBin][nYBin+1];
    
    double values[nXBin][nYBin];
    double w2matrix[nXBin][nYBin]; // sum of squared weights - for error calc
    
    TString name;
    
    


};

#endif // MYH2_HPP