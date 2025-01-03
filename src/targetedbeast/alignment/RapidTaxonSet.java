package targetedbeast.alignment;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;

@Description("A TaxonSet is an ordered set of taxa. The order on the taxa is provided at the time of construction"
		+ " either from a list of taxon objects or an alignment.")
public class RapidTaxonSet extends TaxonSet {

	final public Input<RapidAlignment> alignmentInput = new Input<>("rapidAlignment",
			"alignment where each sequence represents a taxon");

	protected List<String> taxaNames;
	protected List<Taxon> taxonList;

	public RapidTaxonSet() {
	}

	public RapidTaxonSet(final List<Taxon> taxa) {
		taxonsetInput.setValue(taxa, this);
		initAndValidate();
	}

	public RapidTaxonSet(final RapidAlignment alignment) {
		alignmentInput.setValue(alignment, this);
		initAndValidate();
	}

// for testing purposes (Huw)
	public RapidTaxonSet(final String id, final List<Taxon> taxa) {
		setID(id);
		taxonsetInput.setValue(taxa, this);
		initAndValidate();
	}

	@Override
	public void initAndValidate() {

		taxonList = taxonsetInput.get();
		if (alignmentInput.get() != null) {
			if (taxonList.size() > 0) {
				throw new IllegalArgumentException(
						"Only one of taxon and alignment should be specified, not both (id=" + getID() + ").");
			}
			taxaNames = alignmentInput.get().getTaxaNames();
		} else {
			if (taxonList.size() == 0) {
				throw new IllegalArgumentException(
						getID() + ": Either taxon or alignment should be specified (id=" + getID() + ").");
			}
			taxaNames = new ArrayList<>();
			for (final Taxon taxon : taxonList) {
				taxaNames.add(taxon.getID());
			}
		}
	}

	public Set<Taxon> getTaxonSet() {
		final Set<Taxon> unorderedTaxa = new HashSet<>(taxonList);
		return unorderedTaxa;
	}

	/**
	 * @return an unmodifiable list of taxa names as strings.
	 */
	public List<String> asStringList() {
		if (taxaNames == null)
			return null;
		return Collections.unmodifiableList(taxaNames);
	}

	/**
	 * @return the taxa names as a set of strings.
	 */
	public Set<String> getTaxaNames() {
		return new TreeSet<>(taxaNames);
	}

	/**
	 * @return the ID of the i'th taxon.
	 */
	public String getTaxonId(int taxonIndex) {
		return taxaNames.get(taxonIndex);
	}

	/**
	 * return index of given Taxon name
	 * 
	 * @param id
	 * @return -1 if not found
	 */
	public int getTaxonIndex(String id) {
		for (int i = 0; i < taxaNames.size(); i++) {
			if (getTaxonId(i).contentEquals(id))
				return i;
		}
		return -1;
	}

	/**
	 * return Taxon matching the given id
	 * 
	 * @param id
	 * @return null if not found
	 */
	public Taxon getTaxon(String id) {
		for (int i = 0; i < taxonList.size(); i++) {
			Taxon taxon = taxonList.get(i);
			if (taxon.getID().equals(id))
				return taxon;
		}
		return null;
	}

	// convenience methods

	public boolean containsAny(final Collection<String> taxa) {
		final List<String> me = asStringList();
		for (final String taxon : taxa) {
			if (me.contains(taxon)) {
				return true;
			}
		}
		return false;
	}

	public boolean containsAll(final Collection<String> taxa) {
		final List<String> me = asStringList();
		for (final String taxon : taxa) {
			if (!me.contains(taxon)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * @return true if at least 1 member of taxa contained in this set.
	 * @param taxa a collection of taxa
	 */
	public boolean containsAny(final RapidTaxonSet taxa) {
		return containsAny(taxa.asStringList());
	}

	/**
	 * @return true if taxa is a subset of this set
	 * @param taxa
	 */
	public boolean containsAll(final RapidTaxonSet taxa) {
		return containsAll(taxa.asStringList());
	}

	/**
	 * @return number of taxa in this taxon set
	 */
	public int getTaxonCount() {
		if (taxaNames == null)
			return 0;
		return taxaNames.size();
	}

	/**
	 * @return number of taxa in this taxon set
	 * @deprecated Exists only for consistency with method in Alignment. Use
	 *             getTaxonCount() instead.
	 */
	@Deprecated
	public int getNrTaxa() {
		return getTaxonCount();
	}

	@Override
	public String toString() {
		return toString("\t");
	}

	@Override
	public String toString(String indent) {
		final StringBuilder buf = new StringBuilder();
		buf.append(indent).append(getID()).append("\n");
		indent += "\t";
		for (final Taxon taxon : taxonsetInput.get()) {
			buf.append(taxon.toString(indent));
		}
		return buf.toString();
	}
}
