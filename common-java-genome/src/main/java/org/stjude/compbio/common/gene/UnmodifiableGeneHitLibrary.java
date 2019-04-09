package org.stjude.compbio.common.gene;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public class UnmodifiableGeneHitLibrary implements GeneHitLibrary {
	/**
	 * Source library
	 */
	private GeneHitLibrary source;
	
	/**
	 * Constructs an unmodifiable GeneHitLibrary view of another library.
	 */
	public UnmodifiableGeneHitLibrary(GeneHitLibrary source) {
		super();
		this.source = source;
	}

	public Set<String> getAllGeneNames() {
		return Collections.unmodifiableSet(source.getAllGeneNames());
	}

	public Collection<GeneHitGroup> getAllGeneHitGroups() {
		return Collections.unmodifiableCollection(source.getAllGeneHitGroups());
	}

	public GeneHitGroup getGeneHitGroup(String name) {
		return source.getGeneHitGroup(name);
	}

	public List<TxHitGroup> getTxHitGroups(String name) {
		return Collections.unmodifiableList(source.getTxHitGroups(name));
	}

	public TxHitGroup getTxHitGroup(String name) throws IllegalStateException {
		return source.getTxHitGroup(name);
	}

}
