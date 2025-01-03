package targetedbeast.coalescent;

import java.util.Comparator;
import java.util.List;

/**
 * sorts numbers and comparable objects by treating contents of array as a binary beast.tree.
 * KNOWN BUGS: There is a horrible amount of code duplication here!
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 * @version $Id: HeapSort.java,v 1.7 2006/02/20 17:36:23 rambaut Exp $
 */
public class UpdatedSort {


    /**
     * Sorts an array of indices into an array of doubles
     * into increasing order.
     *
     * @param array   an array of doubles
     * @param indices an array of indices to be sorted so that they describe an ascending order of values in array
     */
    public static void sort(double[] array, int[] indices) {

        // ensures we are starting with valid indices
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
        }

        int temp;
        int j, n = indices.length;

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, indices, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = indices[0];
            indices[0] = indices[j];
            indices[j] = temp;
            adjust(array, indices, 1, j);
        }
    }

    // PRIVATE STUFF

    /**
     * helps sort an array of indices pointing into a list of comparable objects.
     * Assumes that indices[lower+1] through to indices[upper] is
     * already in heap form and then puts indices[lower] to
     * indices[upper] in heap form.
     *
     * @param list    a list of comparables
     * @param indices an array of indices pointing into list
     * @param lower   starting index in array to heapify
     * @param upper   end index in array to heapify
     */
    @SuppressWarnings("unchecked")
    private static void adjust(@SuppressWarnings("rawtypes") List<? extends Comparable> list, int[] indices, int lower, int upper) {

        int j, k;
        int temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (list.get(indices[k - 1]).compareTo(list.get(indices[k])) < 0)) {
                k += 1;
            }
            if (list.get(indices[j - 1]).compareTo(list.get(indices[k - 1])) < 0) {
                temp = indices[j - 1];
                indices[j - 1] = indices[k - 1];
                indices[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     *
     * @param array array to sort
     * @param lower lower index of heapify
     * @param upper upper index of heapify
     */
    @SuppressWarnings("unchecked")
    private static void adjust(@SuppressWarnings("rawtypes") Comparable[] array, int lower, int upper) {

        int j, k;
        Comparable<?> temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array[k - 1].compareTo(array[k]) < 0)) {
                k += 1;
            }
            if (array[j - 1].compareTo(array[k - 1]) < 0) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     *
     * @param array      array of objects to sort
     * @param comparator comparator to provide ordering
     * @param lower      lower index of heapify
     * @param upper      upper index of heapify
     */
    @SuppressWarnings("unchecked")
    private static void adjust(Object[] array, @SuppressWarnings("rawtypes") Comparator comparator, int lower, int upper) {

        int j, k;
        Object temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (comparator.compare(array[k - 1], array[k]) < 0)) {
                k += 1;
            }
            if (comparator.compare(array[j - 1], array[k - 1]) < 0) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * helps sort an array of doubles.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     *
     * @param array array of doubles to sort
     * @param lower lower index of heapify
     * @param upper upper index of heapify
     */
    private static void adjust(double[] array, int lower, int upper) {

        int j, k;
        double temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array[k - 1] < array[k])) {
                k += 1;
            }
            if (array[j - 1] < array[k - 1]) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * helps sort an array of doubles.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     *
     * @param array array of doubles to sort ignoring sign
     * @param lower lower index of heapify
     * @param upper upper index of heapify
     */
    private static void adjustAbs(double[] array, int lower, int upper) {

        int j, k;
        double temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (Math.abs(array[k - 1]) < Math.abs(array[k]))) {
                k += 1;
            }
            if (Math.abs(array[j - 1]) < Math.abs(array[k - 1])) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }

    /**
     * helps sort an array of indices into an array of doubles.
     * Assumes that array[lower+1] through to array[upper] is
     * already in heap form and then puts array[lower] to
     * array[upper] in heap form.
     *
     * @param array   array of doubles
     * @param indices array of indices into double array to sort
     * @param lower   lower index of heapify
     * @param upper   upper index of heapify
     */
    private static void adjust(double[] array, int[] indices, int lower, int upper) {

        int j, k;
        int temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array[indices[k - 1]] < array[indices[k]])) {
                k += 1;
            }
            if (array[indices[j - 1]] < array[indices[k - 1]]) {
                temp = indices[j - 1];
                indices[j - 1] = indices[k - 1];
                indices[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }
}

