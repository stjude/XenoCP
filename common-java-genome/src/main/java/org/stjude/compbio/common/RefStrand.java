package org.stjude.compbio.common;

/**
 * Simple implementation of RefStranded that stores ref name and strand as
 * data members.
 * 
 * This is immutable unless extended.
 */
public class RefStrand implements RefStranded {
	private String refName;
	private Strand strand;
	
	/**
	 * Constructor with no initialization--you need to call setters afterwards.
	 */
	protected RefStrand() {
		this(null, null);
	}

	/**
	 * Constructs the AbstractRefStranded with ref name and strand
	 * 
	 * @param refName reference name, can be null, but strongly discouraged
	 * @param refStrand non-null RefStranded
	 * @throws IllegalArgumentException if refStrand is null
	 */
	public RefStrand(String refName, Strand refStrand) {
		this.refName = refName;
		setRefStrand(refStrand);
	}

	/**
	 * Sets the reference name
	 * 
	 * @param refName reference name, can be null, but strongly discouraged
	 */
	protected void setRefName(String refName) {
		this.refName = refName;
	}

	/**
	 * Sets the RefStranded.
	 * 
	 * @param refStrand non-null RefStranded
	 * @throws IllegalArgumentException if refStrand is null
	 */
	protected void setRefStrand(Strand refStrand) {
		if(refStrand == null) throw new IllegalArgumentException(
				"null refStrand");
		this.strand = refStrand;
	}

	public String getRefName() {
		return refName;
	}

	public Strand getStrand() {
		return strand;
	}

	@Override
	public int hashCode() {
		return refName.hashCode() + strand.getOrient();
	}

	@Override
	public boolean equals(Object obj) {
		if(obj == this) return true;
		if(obj == null) return false;
		if(!(obj instanceof RefStranded)) return false;
		
		RefStranded other = (RefStranded)obj;
		if(strand != other.getStrand()) return false;
		if(refName == null) {
			return other.getRefName() == null;
		}
		else {
			return refName.equals(other.getRefName());
		}
	}

	@Override
	public String toString() {
		return refName + "(" + strand + ")";
	}

	
}
