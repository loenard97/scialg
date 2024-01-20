//! Sorting algorithms for onedimensional datasets

/// Returns true only if *arr* is in ascending partial order
///
/// # Example
/// ```
/// use scialg::sort::is_sorted;
///
/// assert!(is_sorted(&vec![1, 2, 3, 4]));
/// assert!(!is_sorted(&vec![3, 1, 2, 4]));
/// assert!(is_sorted(&vec![1, 1, 1, 1]));
/// ```
pub fn is_sorted<T: PartialOrd>(arr: &[T]) -> bool {
    arr.windows(2).all(|w| w[0] <= w[1])
}

/// Sort *arr* inplace using Bubble sort.
///
/// # Example
/// ```
/// use scialg::sort::bubble;
///
/// let mut input = vec![3, 4, 1, 2];
/// let output = vec![1, 2, 3, 4];
///
/// bubble(&mut input);
///
/// assert_eq!(input, output);
/// ```
pub fn bubble<T: PartialOrd>(arr: &mut [T]) {
    let n = arr.len();

    for _ in 0..n - 1 {
        let mut is_sorted = true;
        for i in 0..n - 1 {
            if arr[i] > arr[i + 1] {
                arr.swap(i, i + 1);
                is_sorted = false;
            }
        }
        if is_sorted {
            return;
        }
    }
}

/// Sort *arr* inplace using Insertion sort.
///
/// # Example
/// ```
/// use scialg::sort::insertion;
///
/// let mut input = vec![3, 4, 1, 2];
/// let output = vec![1, 2, 3, 4];
///
/// insertion(&mut input);
///
/// assert_eq!(input, output);
/// ```
pub fn insertion<T: PartialOrd>(arr: &mut [T]) {
    for i in 1..arr.len() {
        let mut j = i;
        while j > 0 && arr[j] < arr[j - 1] {
            arr.swap(j, j - 1);
            j -= 1;
        }
    }
}

/// Sort *arr* inplace using Shell sort.
///
/// # Example
/// ```
/// use scialg::sort::shell;
///
/// let mut input = vec![3, 4, 1, 2];
/// let output = vec![1, 2, 3, 4];
///
/// shell(&mut input);
///
/// assert_eq!(input, output);
/// ```
pub fn shell<T: PartialOrd>(arr: &mut [T]) {
    let mut inc = 1;
    while inc <= arr.len() {
        inc = 3 * inc + 1;
    }

    while inc >= 1 {
        inc /= 3;
        for i in inc..arr.len() {
            let mut j = i;
            while arr[j] < arr[j - inc] {
                arr.swap(j, j - inc);
                j -= inc;
                if j < inc {
                    break;
                }
            }
        }
    }
}
