/****************************************************************
 make-haplotypes.C
 Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/VcfReader.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/Pipe.H"
#include "BOOM/Interval.H"
using namespace std;
using namespace BOOM;

struct VariableRegion {
  String chr;
  Interval interval;
  Vector<Variant> variants;
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  int READLEN;
  void identifyVariableRegions(const Vector<Variant> &,
			       Vector<VariableRegion> &into);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw String("make-haplotypes <in.vcf[.gz]> <genome.2bit> <readlength> <out.fasta>");
  const String vcfFilename=cmd.arg(0);
  const String twoBitFile=cmd.arg(1);
  const int readLen=cmd.arg(2).asInt();
  const String outfile=cmd.arg(3);
  READLEN=readLen;

  // Get list of variable sites
  Vector<Variant> varSites;
  int totalNumSites, numVarSites;
  VcfReader::getVariableSites(vcfFilename,varSites,totalNumSites,numVarSites);
  cerr<<numVarSites<<" variable sites out of "<<totalNumSites<<endl;

  // Identify variable regions
  Vector<VariableRegion> varRegions;
  identifyVariableRegions(varSites,varRegions);
  for(Vector<VariableRegion>::const_iterator cur=varRegions.begin(), end=
	varRegions.end() ; cur!=end ; ++cur) {
    const VariableRegion &region=*cur;
    //cout<<region.interval.length()<<"\t"<<region.variants.size()<<endl;
    /*cout<<region.chr<<"\t"<<region.interval<<"\t";
    for(Vector<Variant>::const_iterator cur=region.variants.begin(), 
	  end=region.variants.end() ; cur!=end ; ++cur)
      cout<<*cur<<"\t";
      cout<<endl;*/
  }

  // Make haplotypes


    // twoBitToFa -seq=chr3 -start=44593 -end=44675 hg19.2bit stdout

  return 0;
}



void Application::identifyVariableRegions(const Vector<Variant> &sites,
					  Vector<VariableRegion> &into)
{
  const int numVarSites=sites.size();
  int currentSite=0;
  while(currentSite<numVarSites) {
    VariableRegion region;
    Variant firstVar=sites[currentSite];
    region.variants.push_back(firstVar);
    region.chr=firstVar.getChr();
    int end=firstVar.getEnd();
    ++currentSite;
    while(currentSite<numVarSites) {
      Variant nextVar=sites[currentSite];
      if(nextVar.getChr()!=region.chr || 
	 nextVar.getBegin()-end>READLEN) break;
      region.variants.push_back(nextVar);
      end=nextVar.getEnd();
      ++currentSite;
    }
    int begin=firstVar.getBegin()-READLEN+1;
    if(begin<0) begin=0;
    end+=READLEN-1; // ### need to check for end of chromosome....?
    region.interval.setBegin(begin);
    region.interval.setEnd(end);
    into.push_back(region);
  }
}



 /*
void Application::identifyVariableRegions(const Vector<Variant> &sites,
					  Vector<VariableRegion> &into)
{
  const int numVarSites=sites.size();
  int currentSite=0;
  while(currentSite<numVarSites) {
    while(currentSite<numVarSites && !sites[currentSite].isIndel()) 
    //while(currentSite<numVarSites && sites[currentSite].isIndel()) 
      ++currentSite;
    if(currentSite>=numVarSites) break;
    Variant firstVar=sites[currentSite];
    ++currentSite;
    while(currentSite<numVarSites && !sites[currentSite].isIndel()) 
    //while(currentSite<numVarSites && sites[currentSite].isIndel()) 
      ++currentSite;
    if(currentSite>=numVarSites) break;
    Variant nextVar=sites[currentSite];
    if(nextVar.getChr()!=firstVar.getChr()) continue;
    int dist=nextVar.getBegin()-firstVar.getEnd();
    cout<<dist<<endl;
    if(dist<10) cout<<firstVar<<"\t"<<nextVar<<endl;
  }
}
 */

