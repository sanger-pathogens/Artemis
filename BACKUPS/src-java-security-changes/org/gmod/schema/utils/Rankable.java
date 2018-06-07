package org.gmod.schema.utils;

/**
 * These domain objects use rank to define an order
 * 
 * @author art
 */
public abstract interface Rankable {

	/**
     * Retrieve the rank, ie order, of this object. Ranks start at zero
     * 
	 * @return the rank
	 */
	public abstract int getRank();
}
