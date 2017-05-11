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
#include "BOOM/Map.H"
#include "BOOM/Time.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 VariableRegion
 ****************************************************************/
struct VariableRegion {
  String chr;
  Interval interval;
  Vector<Variant> variants; // indels only
  bool contains(const Variant &);
};

/****************************************************************
 HaplotypeRecord
 ****************************************************************/
struct HaplotypeRecord {
  Vector<int> alleles;
  String hash();
};

/****************************************************************
 Application
 ****************************************************************/
class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  int READLEN;
  Map<String,int> chromLen;
  bool wantPrependChr;
  void identifyVariableRegions(const Vector<Variant> &,
			       Vector<VariableRegion> &into);
  int numIndels(const Vector<Variant> &);
  void getChromLengths(const String &twoBitFile,Map<String,int> &into);
  void prependChr(Vector<Variant> &);
  void makeNgenome(Vector<Variant> &,const String &twoBitFile,
		   const String &outfile);
  void partitionByChr(const Vector<Variant> &,
		      Map<String,Vector<Variant> > &into);
  void getSeq(const String &twoBitFile,const String &chr,
	      int begin,int end,String &into);
  void getSeq(const String &twoBitFile,const VariableRegion &region,
	      String &into);
  void setNs(Vector<Variant> &,String &seq);
  void getIndels(const Vector<Variant> &from,Vector<Variant> &into);
  void makeHaplotypes(const Vector<VariableRegion> &,const String &vcfFile,
		      const String &twoBitFile,const String &outfile);
  void makeHaplotypes(const VariableRegion &region,
		      const String &ref,
		      VcfReader &,
		      const String &twoBitFile,
		      const String &outfile);
  void makeHapRecs(const Vector<VariantAndGenotypes> &,
		   Vector<HaplotypeRecord> &into);
  //void uniquify(Vector<HaplotypeRecord> &);
  void personalize(String &seq,const HaplotypeRecord &,
		   const Vector<VariantAndGenotypes> &,int regionStart);
};


/****************************************************************
 main()
 ****************************************************************/
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



bool VariableRegion::contains(const Variant &v)
{
  return chr==v.getChr() && interval.contains(v.getPos());
}



String HaplotypeRecord::hash()
{
  String h;
  for(Vector<int>::const_iterator cur=alleles.begin(), end=alleles.end() ;
      cur!=end ; ++cur)
    h+=String(" ")+(*cur);
  return h;
}



Application::Application()
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"c");
  if(cmd.numArgs()!=4)
    throw String("\n\
make-haplotypes [opt] <in.vcf[.gz]> <genome.2bit> <readlength> <out.fasta>\n\
  NOTE: twoBitToFa and twoBitInfo must be on your $PATH\n\
  NOTE: If using bowtie, use the --np 0 and --n-ceil L,0,0.5 options\n\
        to not penalize N's\n\
  -c = prepend \"chr\" to chrom names in VCF\n\
");
  const String vcfFilename=cmd.arg(0);
  const String twoBitFile=cmd.arg(1);
  const int readLen=cmd.arg(2).asInt();
  const String outfile=cmd.arg(3);
  READLEN=readLen;
  wantPrependChr=cmd.option('c');

  // Get chromosome lengths
  getChromLengths(twoBitFile,chromLen);

  // Make sure outfile doesn't exist
  if(File::exists(outfile)) system((String("rm ")+outfile).c_str());

  // Get list of variable sites
  cout<<getDateAndTime()<<endl;
  cout<<"Getting list of variable sites"<<endl;
  Vector<Variant> varSites;
  int totalNumSites, numVarSites;
  VcfReader::getVariableSites(vcfFilename,varSites,totalNumSites,numVarSites);
  if(wantPrependChr) prependChr(varSites);
  cout<<numVarSites<<" variable sites out of "<<totalNumSites<<endl;

  // Identify variable regions containing indels
  cout<<getDateAndTime()<<endl;
  cout<<"Identifying variable regions"<<endl;
  Vector<VariableRegion> varRegions;
  identifyVariableRegions(varSites,varRegions);
  for(Vector<VariableRegion>::const_iterator cur=varRegions.begin(), end=
	varRegions.end() ; cur!=end ; ++cur) {
    const VariableRegion &region=*cur;
  }

  // Make haplotypes
  cout<<getDateAndTime()<<endl;
  cout<<"Making indel haplotypes"<<endl;
  makeHaplotypes(varRegions,vcfFilename,twoBitFile,outfile);

  // Make N-genome (genomic sequence with N's substituted for SNPs)
  cout<<getDateAndTime()<<endl;
  cout<<"Making N-genome"<<endl;
  makeNgenome(varSites,twoBitFile,outfile);

  cout<<getDateAndTime()<<endl;
  cout<<"[done]"<<endl;
  return 0;
}



void Application::identifyVariableRegions(const Vector<Variant> &allSites,
					  Vector<VariableRegion> &into)
{
  Vector<Variant> indels;
  getIndels(allSites,indels);
  const int numIndels=indels.size();
  int currentSite=0;
  while(currentSite<numIndels) {
    VariableRegion region;
    Variant firstVar=indels[currentSite];
    region.variants.push_back(firstVar);
    region.chr=firstVar.getChr();
    int end=firstVar.getEnd();
    ++currentSite;
    while(currentSite<numIndels) {
      Variant nextVar=indels[currentSite];
      if(nextVar.getChr()!=region.chr || 
	 nextVar.getBegin()-end>=READLEN) break;
      region.variants.push_back(nextVar);
      end=nextVar.getEnd();
      ++currentSite;
    }
    int begin=firstVar.getBegin()-READLEN+1;
    if(begin<0) begin=0;
    end+=READLEN-1;
    if(!chromLen.isDefined(region.chr))
      throw String("chromosome ")+region.chr+" is not present in 2bit file";
    const int L=chromLen[region.chr];
    if(end>L) end=L;
    region.interval.setBegin(begin);
    region.interval.setEnd(end);
    into.push_back(region);
  }
}



int Application::numIndels(const Vector<Variant> &variants)
{
  int n=0;
  for(Vector<Variant>::const_iterator cur=variants.begin(), end=
	variants.end() ; cur!=end ; ++cur)
    if((*cur).isIndel()) ++n;
  return n;
}



void Application::getChromLengths(const String &twoBitFile,
				  Map<String,int> &into)
{
  String cmd=String("twoBitInfo ")+twoBitFile+" stdout";
  Pipe pipe(cmd,"r");
  while(!pipe.eof()) {
    String line=pipe.getline(); line.trimWhitespace();
    Vector<String> fields;
    line.getFields(fields);
    if(fields.size()!=2) continue;
    String chr=fields[0];
    const int L=int(fields[1]);
    into[chr]=L;
  }
  pipe.close();
}



void Application::prependChr(Vector<Variant> &variants)
{
  for(Vector<Variant>::iterator cur=variants.begin(), end=variants.end() ;
      cur!=end ; ++cur) {
    Variant &variant=*cur;
    String &chr=variant.getChr();
    chr=String("chr")+chr;
  }
}



void Application::partitionByChr(const Vector<Variant> &variants,
				 Map<String,Vector<Variant> > &into)
{
  for(Vector<Variant>::const_iterator cur=variants.begin(),
	end=variants.end() ; cur!=end ; ++cur) {
    const Variant &v=*cur;
    const String &chr=v.getChr();
    into[chr].push_back(v);
  }
}



void Application::getSeq(const String &twoBitFile,const String &chr,
			 int begin,int end,String &into)
{
  String cmd=String("twoBitToFa -seq=")+chr+" -start="+begin+" -end="+end+
    " "+twoBitFile+" stdout";
  Pipe pipe(cmd,"r");
  String defline=pipe.getline();
  while(!pipe.eof()) into+=pipe.getline();
}



void Application::getSeq(const String &twoBitFile,
			 const VariableRegion &region,
			 String &into)
{
  getSeq(twoBitFile,region.chr,region.interval.getBegin(),
	 region.interval.getEnd(),into);
}



void Application::makeNgenome(Vector<Variant> &variants,
			      const String &twoBitFile,
			      const String &outfile)
{
  // Partition SNPs by chr
  Map<String,Vector<Variant> > byChr;
  partitionByChr(variants,byChr);

  // Process twoBitFile by chr
  FastaWriter writer;
  Set<String> chroms;
  chromLen.getKeys(chroms);
  for(Set<String>::const_iterator cur=chroms.begin(), end=chroms.end() ;
      cur!=end ; ++cur) {
    const String &chr=*cur;
    const int L=chromLen[chr];
    String seq;
    getSeq(twoBitFile,chr,0,L,seq);
    if(byChr.isDefined(chr)) setNs(byChr[chr],seq);
    writer.appendToFasta(String(">")+chr,seq,outfile);
  }
}



void Application::setNs(Vector<Variant> &variants,String &seq)
{
  for(Vector<Variant>::const_iterator cur=variants.begin(),
	end=variants.end() ; cur!=end ; ++cur) {
    Variant v=*cur;
    if(v.isIndel()) continue;
    const int pos=v.getPos()-1;
    if(v.getAllele(0)[0]!=toupper(seq[pos]))
      cout<<"ref mismatch: "+seq.substring(pos-1,3)+" vs "+v.getAllele(0)
	  <<endl;
    seq[pos]='N';
  }
}



void Application::getIndels(const Vector<Variant> &from,
			    Vector<Variant> &into)
{
  for(Vector<Variant>::const_iterator cur=from.begin(), end=from.end() ;
      cur!=end ; ++cur)
    if((*cur).isIndel()) into.push_back(*cur);
}



void Application::makeHaplotypes(const Vector<VariableRegion> &regions,
				 const String &vcfFile,
				 const String &twoBitFile,
				 const String &outfile)
{
  VcfReader reader(vcfFile);
  for(Vector<VariableRegion>::const_iterator cur=regions.begin(),
	end=regions.end() ; cur!=end ; ++cur) {
    const VariableRegion &region=*cur;
    String ref;
    getSeq(twoBitFile,region,ref);
    makeHaplotypes(region,ref,reader,twoBitFile,outfile);
  }
}



void Application::makeHapRecs(const Vector<VariantAndGenotypes> &recs,
			      Vector<HaplotypeRecord> &hapRecs)
{
  if(recs.size()==0) INTERNAL_ERROR;
  const VariantAndGenotypes &firstRec=recs[0];
  const int numGenotypes=firstRec.genotypes.size();
  for(int g=0 ; g<numGenotypes ; ++g) {
    const int numAlleles=firstRec.genotypes[g].numAlleles();
    for(int a=0 ; a<numAlleles ; ++a) {
      HaplotypeRecord haprec;
      for(Vector<VariantAndGenotypes>::const_iterator cur=recs.begin(),
	    end=recs.end() ; cur!=end ; ++cur)
	haprec.alleles.push_back((*cur).genotypes[g][a]);
      hapRecs.push_back(haprec);
    }
  }
}



/*
void Application::uniquify(Vector<HaplotypeRecord> &hapRecs)
{
  Vector<HaplotypeRecord> temp;
  Set<String> seen;
  for(Vector<HaplotypeRecord>::const_iterator cur=hapRecs.begin(),
	end=hapRecs.end() ; cur!=end ; ++cur) {
    const HaplotypeRecord &rec=*cur;
    const String hash=rec.hash();
    if(seen.isMember(hash)) continue;
    seen+=hash;
    temp.push_back(rec);
  }
  hapRecs=temp;
}
*/


void Application::makeHaplotypes(const VariableRegion &region,
				 const String &ref,
				 VcfReader &reader,
				 const String &twoBitFile,
				 const String &outfile)
{
  Vector<VariantAndGenotypes> records;
  VariantAndGenotypes rec;

  // Advance to the first variant in this region
  while(reader.currentVariant(rec)) {
    if(wantPrependChr) 
      rec.variant.getChr()=String("chr")+rec.variant.getChr();
    if(!region.contains(rec.variant)) {
      reader.advance();
      continue; }
    else break;
  }

  // Collect variants in this region
  while(reader.currentVariant(rec)) {
    if(wantPrependChr) 
      rec.variant.getChr()=String("chr")+rec.variant.getChr();
    if(region.contains(rec.variant)) {
      records.push_back(rec);
      reader.advance(); }
    else break;
  }

  // Process region to produce haplotypes
  if(records.empty()) throw "No records found for region";
  Vector<HaplotypeRecord> hapRecs;
  makeHapRecs(records,hapRecs);
  FastaWriter writer;
  Set<String> seen;
  String baseDef=String(">hap:")+region.chr+":"+region.interval.getBegin()
    +"-"+region.interval.getEnd();
  int hapID=1;
  for(Vector<HaplotypeRecord>::const_iterator cur=hapRecs.begin(), 
	end=hapRecs.end() ; cur!=end ; ++cur) {
    const HaplotypeRecord &hapRec=*cur;
    String alt=ref;
    personalize(alt,hapRec,records,region.interval.getBegin());
    if(seen.isMember(alt)) continue;
    seen+=alt;
    String def=baseDef+"_"+hapID;
    writer.appendToFasta(def,alt,outfile);
    ++hapID;
  }
}



void Application::personalize(String &seq,const HaplotypeRecord &hapRec,
			      const Vector<VariantAndGenotypes> &variants,
			      int regionStart)
{
  //TRACE
  int numVariants=variants.size();
  //TRACE
  if(hapRec.alleles.size()!=numVariants) INTERNAL_ERROR;
  //TRACE
  for(int i=numVariants-1 ; i>=0 ; --i) {
    //TRACE
    const Variant &variant=variants[i].variant;
    //TRACE
    if(variant.isIndel()) {
      //TRACE
      const int allele=hapRec.alleles[i];
      //TRACE
      if(allele==0) continue; // reference
      //TRACE
      const int refLen=variant.getAllele(0).length();
      //TRACE
      const String &altAllele=variant.getAllele(allele);
      //TRACE
      //cout<<"indel "<<variant.getPos()<<" "<<regionStart<<" "<<variant.getPos()-regionStart<<" "<<refLen<<endl;
      seq.replaceSubstring(variant.getPos()-1-regionStart,refLen,altAllele);
      //TRACE
    }
    else {
      //cout<<"SNP "<<variant.getPos()<<" "<<regionStart<<" "<<variant.getPos()-regionStart<<endl;
      seq[variant.getPos()-1-regionStart]='N';
    }
    //TRACE
  }
  //TRACE
}



