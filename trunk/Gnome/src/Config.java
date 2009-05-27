/**
 * 
 */

/**
 * @author xzhou
 *
 */

/* *
 * casectrl_77SNP.fasta 77 3000 true 1000000000 50 2120 8 5 5 -3 9 freq
 * */
public class Config {
	public String caseControlFileName;
	public int nSNPs;
	public int nRecords;
	public Boolean allGoodSNP;
	public int rsquarePrecision;
	public int pValuePrecision;
	public int nIndividualInCase;
	public int sigSNPIndex1;
	public int sigSNPIndex2;
	public int approachChoice;
	
	public Config(String confFileName)
	{
		readConfigFile(confFileName);
	}
	
	public void readConfigFile(String confFileName)
	{
		
	}
}
