import java.io.Serializable;
import java.util.Enumeration;
import java.util.Vector;

import jp.ac.kobe_u.cs.cream.IntVariable;

enum Allele implements Serializable {
    A,T,G,C,U  //U means unknown
}

class Entry implements Serializable {
	private long iCaseA;
	public long getCaseA() {
		return iCaseA;
	}
	private long iCaseT;
	public long getCaseT() {
		return iCaseT;
	}
	
	public Entry(long iCaseA,long iCaseT) {
		this.iCaseA = iCaseA;
		this.iCaseT = iCaseT;
	}

	public boolean equal(Entry entry) {
		if (this.iCaseA == entry.iCaseA && this.iCaseT == entry.iCaseT) 
			return true;
		else
			return false;
	}
	
	public String toString() {
		return "iCaseA="+iCaseA+" iCaseT="+iCaseT;
	}
}

class Entries implements Serializable {
	public Vector<Entry> entries = new Vector<Entry>();
	
	public void appendWithoutDuplicate(Entry entry) {
		for (Enumeration<Entry> e = entries.elements();e.hasMoreElements();) {
			Entry en = e.nextElement();
			if (en.equal(entry)) return;
		}
		entries.add(entry);
	}
	
	public int getEntryNum() {
		return entries.size();
	}
}

class ConstraintForRSquare implements Serializable { //Constraint for rsquare
	public IntVariable vAC;
	public IntVariable vAD;
	public IntVariable vBC;
	public IntVariable vBD;
	
	public int iAC,iAD,iBC,iBD;
	
	public ConstraintForRSquare(IntVariable vAC,IntVariable vAD,IntVariable vBC,IntVariable vBD) {
		this.vAC = vAC;
		this.vAD = vAD;
		this.vBC = vBC;
		this.vBD = vBD;
	}
}

class Segment implements Serializable { 
	public Segment(int begin, int length) {
		this.begin = begin;
		this.length = length;
	}
	public String toString() {
		return ""+this.begin+"-"+this.length;
	}
	int begin;
	int length;
}

class Record implements Serializable {
	public Vector<Allele> record = new Vector<Allele>();
	public int count; //number of record
	public int min; //min of count
	public int max; //max of count
	public Vector<Segment> segments = new Vector<Segment>();
	public ConstraintForRecord constraint;
	private int begin;
	private boolean free = false;
	private int lb = Integer.MIN_VALUE;
	public int getLB() {
		return lb;
	}

	private int ub = Integer.MAX_VALUE;
	private boolean removable = false;
	
	public int getUB() {
		return ub;
	}

	public boolean isFree() {
		return free;
	}

	public Record(int count,Allele[] r,int begin) {
		this.count = count;
		
		for (int i=0;i<r.length;i++) {
			this.record.add(r[i]);
		}
		
		this.begin = begin;
	}
	
	public int getSNPIDByIndex(int index) {
		int i = 0;
		for (Enumeration<Segment> segs = segments.elements();segs.hasMoreElements();) {
			Segment seg = segs.nextElement();
			if (i+seg.length <= index) {
				i += seg.length;
			} else {
				return seg.begin + (index - i);
			}
		}
		return -1;
	}

	public Record(int count, Vector<Allele> alleles, int begin) {
		this.count = count;
		this.record.addAll(alleles);
		this.begin = begin;
	}
	
	public String toString(Allele[] refA,Allele[] refT) {
		int i = 0;
		String minmax = "";
		if (this.min != this.max) minmax = "("+this.min+","+this.max+")";
		String ret = ""+this.count+minmax+"   [";
		for (Enumeration<Allele> a = record.elements();a.hasMoreElements();) {
			Allele ele = a.nextElement();
			if (ele == Allele.A) 
				ret += refA[i]+" ";
			else
				ret += refT[i]+" ";
			i++;
		}
		ret += "]";
		return ret;
	}
	
	public String toString() {
		int i = 0;
		String minmax = "";
		if (this.min != this.max) minmax = "("+this.min+","+this.max+")";
		String ret = ""+this.count+minmax+"   [";
		
		for (Enumeration<Allele> a = record.elements();a.hasMoreElements();) {
			Allele ele = a.nextElement(); 
			ret += ele+" ";
			i++;
		}
		ret += "]";
		return ret;
	}

	public int getBeginIndex() {
		return begin;
	}

	public int getIndexBySNPID(int id) {
		int i = 0;
		for (Enumeration<Segment> segs = segments.elements();segs.hasMoreElements();) {
			Segment seg = segs.nextElement();
			if (id >= seg.begin && id < seg.begin+seg.length) {
				return i+(id-seg.begin);
			}
			i += seg.length;
		}
		return -1;
	}

	public void setFree(boolean b) {
		this.free  = b;
	}

	public void setLB(int lb) throws Exception {
		this.lb = Math.max(this.lb, lb);
		if (this.lb < 0) throw new Exception("negative lower bound "+lb);
	}

	public void setUB(int ub) throws Exception {
		this.ub = Math.min(this.ub, ub);
		if (this.ub < 0) throw new Exception("negative upper bound "+ub);
	}

	public boolean isRemovable() {
		return this.removable ;
	}

	public void setRemovable(boolean b) {
		this.removable = b;
	}

	public String toStringWithDetails() {
		String ret = "";
		for (Enumeration<Allele> a = record.elements();a.hasMoreElements();) {
			Allele ele = a.nextElement(); 
			ret += ele+" ";
		}
		ret += "\n";
		for (Enumeration<Segment> a = segments.elements();a.hasMoreElements();) {
			Segment ele = a.nextElement(); 
			for (int i=ele.begin;i<ele.begin+ele.length;i++) {
				ret += "" + i + " ";
			}
		}
		return ret;
	}
}

class ConstraintForRecord implements Serializable { //Constraint for record: the constraint between count and frequency
	public IntVariable count;
	public int value;
	
	public ConstraintForRecord(IntVariable count) {
		this.count = count;
		count.ge(0);
	}
}

class SNP implements Serializable {
	private int ID;
	private String name;
	private int iCASE;

	public Vector<RSquare> rsquares = new Vector<RSquare>();
	private double pvalue;
	public double getPvalue() {
		return pvalue;
	}

	private int iCaseAFreq=0, iCaseTFreq=0, iControlAFreq=0, iControlTFreq=0;
	public int getControlTFreq() {
		return iControlTFreq;
	}

	public int getControlAFreq() {
		return iControlAFreq;
	}

	public int getCaseTFreq() {
		return iCaseTFreq;
	}

	public int getCaseAFreq() {
		return iCaseAFreq;
	}
	
	public int getCasePopulation() {
		return this.iCASE;
	}
	
	public int getControlPopulation() {
		return this.alleles.size()-this.iCASE;
	}
	
	public int getPopulation() {
		return this.alleles.size();
	}
	
	public int getID() {
		return ID;
	}

	public String getName() {
		return name;
	}

	private Vector<Allele> alleles = new Vector<Allele>();
	private Allele MAF;
	private boolean badSNP = false;
	public Entries entries;
	private boolean involved = false;
	private boolean sigSNP = false;
	
	public boolean isSigSNP() {
		return sigSNP;
	}

	public boolean isInvolved() {
		return involved;
	}

	public boolean isBadSNP() {
		return badSNP;
	}

	public Allele getMAF() {
		return MAF;
	}

	public Vector<Allele> getAlleles() {
		return alleles;
	}

	public SNP(String name,Vector<Allele> alleles,int iCASE) {
		this.name = name;
		this.alleles.addAll(alleles);
		this.iCASE = iCASE;
		
		updateStatus();
	}
	
	//update the values of all variables
	private void updateStatus() {
		int i = 0;
		for (Enumeration<Allele> a = alleles.elements();a.hasMoreElements();) {
			Allele allele = a.nextElement();
			if (i<iCASE) {
				if (allele == Allele.A) this.iCaseAFreq++;
				else if (allele == Allele.T) this.iCaseTFreq++;
			} else {
				if (allele == Allele.A) this.iControlAFreq++;
				else if (allele == Allele.T) this.iControlTFreq++;
			}
			i++;
		}
		
		if (iCaseAFreq >= iCaseTFreq) {
			MAF = Allele.T;
//			if (iControlAFreq >= iControlTFreq)
//				MAF = Allele.T;
//			else
//				MAF = Allele.U;
		} else {
			MAF = Allele.A;
//			if (iControlAFreq < iControlTFreq)
//				MAF = Allele.A;
//			else
//				MAF = Allele.U;
		}
		
		this.pvalue = computePVALUE();
	}

	public SNP(String name,int ID,Allele[] alleles,int iCASE) {
		this.name = name;
		this.ID = ID;
		for (int i=0;i<alleles.length;i++)  this.alleles.add(alleles[i]);
		this.iCASE = iCASE;
		
		updateStatus();
	}
	
	//calculate RSQAURE
	private double computePVALUE() {
		int iAC=0,iAD=0,iBC=0,iBD=0;
		int count = 0;
		for (Enumeration<Allele> a = alleles.elements();a.hasMoreElements();) {
			Allele allele = a.nextElement();
			if (count < this.iCASE) {
				if (allele == Allele.A) iAC++;
				else if (allele == Allele.T) iBC++;
			} else {
				if (allele == Allele.A) iAD++;
				else if (allele == Allele.T) iBD++;
			}
			count++;
		}
		
		return calculateChiSquare(iAC,iAD,iBC,iBD);
//		int total = iAC + iAD + iBC + iBD;
//		long x11,x12,x21,x22,p1,q1,p2,q2;
//		x11 = iAC;
//		x12 = iAD;
//		x21 = iBC;
//		x22 = iBD;
//		
//		p1 = x11 + x12;
//		q1 = x11 + x21;
//		p2 = x21 + x22;
//		q2 = x12 + x22;
//		
//		//System.out.println(""+Math.pow(x11*total-p1*q1,2));
//		//System.out.println(""+(p1*q1*p2*q2));
//		double pvalue = 0;
//		long div = p1*q1*p2*q2;
//		if (div == 0) 
//			pvalue = -1;
//		else
//			pvalue = Math.pow(x11*total-p1*q1,2) / (p1*q1*p2*q2);
//		
//		return pvalue;
	}
	
	private double calculateChiSquare(int iAC, int iAD, int iBC, int iBD) {
		if (iAC < 0 || iAD < 0 || iBC < 0 || iBD < 0) {
			System.err.println("error calculateRSQUARE("+iAC+","+iAD+","+iBC+","+iBD+")");
			System.exit(0);
		}
		
		long total = iAC + iAD + iBC + iBD;
		double x11,x12,x21,x22,p1,q1,p2,q2,chisquare;
		x11 = iAC*1.0;
		x12 = iAD*1.0;
		x21 = iBC*1.0;
		x22 = iBD*1.0;
		
		p1 = x11 + x12;
		q1 = x11 + x21;
		p2 = x21 + x22;
		q2 = x12 + x22;
		double div = p1*p2*q1*q2;
		if (div == 0) 
			chisquare = -1;
		else
			//chisquare = Math.pow((x11*x22-x12*x21), 2)*(p1+p2)/div;
			chisquare = Math.pow((x11*x22-x12*x21), 2)/div;
		
		return chisquare;
	}

	public void appendRSquare(RSquare rsquare) {
		rsquares.add(rsquare);
	}

	public void setBadSNP(boolean b) {
		badSNP  = b;
	}

	public void setInvolved(boolean involved) {
		this.involved  = involved;
	}

	public void setSigSNP(boolean b) {
		this.sigSNP  = b;
	}

	public void setPvalue(double p) {
		this.pvalue = p;
	}

	public void setCaseAFreq(int c) {
		this.iCaseAFreq = c;
		this.iCASE = this.iCaseAFreq+this.iCaseTFreq;
	}

	public void setCaseTFreq(int c) {
		this.iCaseTFreq = c;
		this.iCASE = this.iCaseAFreq+this.iCaseTFreq;
	}

	public Double getpA() {
		if (this.iCaseAFreq >= this.iCaseTFreq)
			return this.iCaseAFreq*1.0/this.iCASE;
		else
			return this.iCaseTFreq*1.0/this.iCASE;
	}
}

class Block implements Serializable {
	private int beginIndex;
	public int getBeginIndex() {
		return beginIndex;
	}

	private int endIndex;

	public int getEndIndex() {
		return endIndex;
	}

	public Block(SNP[] snps, int begin, int length) {
		for (int i=begin;i<begin+length;i++) {
			this.snps.add(snps[i]);
		}
		this.beginIndex = begin;
		this.endIndex = begin+length-1;
	}

	public Vector<SNP> snps = new Vector<SNP>();
	public void generatePairWiseBlock(SNP[] snps, int index1, int index2) {
		this.snps.clear();
		this.snps.add(snps[index1]);
		this.snps.add(snps[index2]);
	}
}

public class RSquare implements Serializable {
	public SNP snp1;
	public SNP snp2;
	
	public int iAC=0, iAD=0, iBC=0, iBD=0;
	public double rsquare; //r2
	public double r; //r
	public ConstraintForRSquare constraint;
	private boolean badRSquare = false;
	private boolean recovered = false;
	private boolean signRecovered = false;
	public double realr2;
	public double pAB;
	public double pA;
	public double pB;
	public char sign;
	public double v1; //r candidate 1
	public double v2; //r candidate 2
	
	public boolean isSignRecovered() {
		return signRecovered;
	}

	public boolean isRecovered() {
		return recovered;
	}

	public boolean isBadRSquare() {
		return badRSquare;
	}

	public int getTotal() {
		return iAC+iAD+iBC+iBD;
	}

	public RSquare(SNP snp1, SNP snp2) {
		// calculate rsquare
		
		Vector<Allele> alleles1 = snp1.getAlleles();
		Vector<Allele> alleles2 = snp2.getAlleles();
		int iCASE = snp1.getCasePopulation();
		//int iCASE = snp1.getPopulation();
		for (Enumeration<Allele> e1 = alleles1.elements(),e2 = alleles2.elements(); 
			e1.hasMoreElements();) {
			iCASE--;
			if (iCASE == -1) break;
			Allele a1 = e1.nextElement();
			Allele a2 = e2.nextElement();
			if (a1 == Allele.A && a2 == Allele.A) iAC++;
			else if (a1 == Allele.A && a2 == Allele.T) iAD++;
			else if (a1 == Allele.T && a2 == Allele.A) iBC++;
			else if (a1 == Allele.T && a2 == Allele.T) iBD++;
		}
		
		rsquare = calculateRSQUARE(iAC,iAD,iBC,iBD);
		r = iAC*iBD - iAD*iBC;
		//r = r/Math.abs(r);
		int s1 = 1, s2 = 1;
		if (snp1.getCaseAFreq() > snp1.getCaseTFreq()) s1 = -1;
		if (snp2.getCaseAFreq() > snp2.getCaseTFreq()) s2 = -1;
		if (r>0) r = Math.sqrt(rsquare)*s1*s2;
		else if (r<0) r = -Math.sqrt(rsquare)*s1*s2;
		//System.out.println("r="+r);
		this.snp1 = snp1;
		this.snp2 = snp2;
		//System.out.println("rsquare="+rsquare+" for snp_"+snp1.getID()+"_"+snp2.getID());
	}
	
	//calculate RSQAURE
	public static double calculateRSQUARE(long iAC, long iAD, long iBC, long iBD) {
		if (iAC < 0 || iAD < 0 || iBC < 0 || iBD < 0) {
			System.err.println("error calculateRSQUARE("+iAC+","+iAD+","+iBC+","+iBD+")");
			System.exit(0);
		}
		
		long total = iAC + iAD + iBC + iBD;
		long x11,x12,x21,x22,p1,q1,p2,q2;
		x11 = iAC;
		x12 = iAD;
		x21 = iBC;
		x22 = iBD;
		
		p1 = x11 + x12;
		q1 = x11 + x21;
		p2 = x21 + x22;
		q2 = x12 + x22;
		
		//System.out.println(""+Math.pow(x11*total-p1*q1,2));
		//System.out.println(""+(p1*q1*p2*q2));
		double rsquare = 0;
		long div = p1*q1*p2*q2;
		if (div == 0) 
			rsquare = -1;
		else
			rsquare = Math.pow(x11*total-p1*q1,2) / (p1*q1*p2*q2);
		
		return rsquare;
	}

	public void setBadRSquare(boolean b) {
		badRSquare  = b;
	}

	public void setRecovered(boolean b) {
		this.recovered  = b;
	}

	public void setSignRecovered(boolean b) {
		this.signRecovered  = b;
		//System.out.println("signRecovered,r="+this.r+",r2="+this.rsquare);
	}

	public double getpAB() {
		if (snp1.getCaseAFreq()>=snp1.getCaseTFreq()) {
			if (snp2.getCaseAFreq()>=snp2.getCaseTFreq())
				return iAC*1.0/snp1.getCasePopulation();
			else 
				return iAD*1.0/snp1.getCasePopulation();
		} else {
			if (snp2.getCaseAFreq()>=snp2.getCaseTFreq())
				return iBC*1.0/snp1.getCasePopulation();
			else 
				return iBD*1.0/snp1.getCasePopulation();
		}
	}

	public double getSign() {
//		double s1 = 1,s2 = 1;
//		if (rs.snp1.getCaseAFreq()<rs.snp1.getCaseTFreq()) {
//			s1 = -1;
//		}
//		if (rs.snp2.getCaseAFreq()<rs.snp2.getCaseTFreq()) {
//			s2 = -1;
//		}
//		double sign = s1*s2*rs.r;
//		if (sign != 0) sign = sign/Math.abs(sign);
//		return sign;
//		int D = iAC*snp1.getCasePopulation()-snp1.getCaseAFreq()*snp2.getCaseAFreq();
//		return D;
		
		int D = 0;
		if (snp1.getCaseAFreq()>=snp1.getCaseTFreq()) {
			if (snp2.getCaseAFreq()>=snp2.getCaseTFreq())
				D = iAC*snp1.getCasePopulation()-snp1.getCaseAFreq()*snp2.getCaseAFreq();
			else 
				D = iAD*snp1.getCasePopulation()-snp1.getCaseAFreq()*snp2.getCaseTFreq();
		} else {
			if (snp2.getCaseAFreq()>=snp2.getCaseTFreq())
				D = iBC*snp1.getCasePopulation()-snp1.getCaseTFreq()*snp2.getCaseAFreq();
			else 
				D = iBD*snp1.getCasePopulation()-snp1.getCaseTFreq()*snp2.getCaseTFreq();
		}
		if (D != 0) D = D/Math.abs(D);
		return D;
	}

	public double getpB() {
//		if (snp2.getCaseAFreq()>=snp2.getCaseTFreq()) {
//			return snp2.getCaseAFreq()*1.0/snp2.getCasePopulation();
//		} else {
//			return snp2.getCaseTFreq()*1.0/snp2.getCasePopulation();
//		}
		return this.pB;
	}
	
	public double getpA() {
//		if (snp1.getCaseAFreq()>=snp1.getCaseTFreq()) {
//			return snp1.getCaseAFreq()*1.0/snp1.getCasePopulation();
//		} else {
//			return snp1.getCaseTFreq()*1.0/snp1.getCasePopulation();
//		}
		return this.pA;
	}

	public void setSignRecovered(boolean sign, boolean b) {
		if (sign) {
			this.pAB = v1;
		} else {
			this.pAB = v2;
		}
		this.signRecovered = b;
		
		if (b) {
			if ((sign && this.r <= 0) || (!sign && this.r >= 0)) {
				System.err.println("error detected: incorrect sign");
			}
		}
		System.out.println("pair ("+this.snp1.getID()+","+this.snp2.getID()+") recovered,r2="+this.rsquare);
	}

	public double getRecoveredSign() throws Exception {
		if (!this.signRecovered) {
			throw new Exception("sign is not recovered");
		}
		if (this.pAB == v1) return 1;
		else if (this.pAB == v2) return -1;
		else throw new Exception("invalid sign");
	}
}