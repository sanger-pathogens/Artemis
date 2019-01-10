package org.gmod.schema.utils;

public class PeptideProperties {
	
	private String mass;
	
	private String aminoAcids;
	
	private String isoelectricPoint;
	
	private String charge;
	
	private String signalPeptide;
	
	private String transmembraneDomain;
	
	private String gpiAnchor;

	public String getAminoAcids() {
		return aminoAcids;
	}

	public void setAminoAcids(String aminoAcids) {
		this.aminoAcids = aminoAcids;
	}

	public String getCharge() {
		return charge;
	}

	public void setCharge(String charge) {
		this.charge = charge;
	}

	public String getGpiAnchor() {
		return gpiAnchor;
	}

	public void setGpiAnchor(String gpiAnchor) {
		this.gpiAnchor = gpiAnchor;
	}

	public String getIsoelectricPoint() {
		return isoelectricPoint;
	}

	public void setIsoelectricPoint(String isoelectricPoint) {
		this.isoelectricPoint = isoelectricPoint;
	}

	public String getMass() {
		return mass;
	}

	public void setMass(String mass) {
		this.mass = mass;
	}

	public String getSignalPeptide() {
		return signalPeptide;
	}

	public void setSignalPeptide(String signalPeptide) {
		this.signalPeptide = signalPeptide;
	}

	public String getTransmembraneDomain() {
		return transmembraneDomain;
	}

	public void setTransmembraneDomain(String transmembraneDomain) {
		this.transmembraneDomain = transmembraneDomain;
	}
	
	
}
