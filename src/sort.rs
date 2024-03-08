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

/// Sort *arr* inplace using exchange sort.
///
/// # Example
/// ```
/// use scialg::sort::exchange;
///
/// let mut input = vec![3, 4, 1, 2];
/// let output = vec![1, 2, 3, 4];
///
/// exchange(&mut input);
///
/// assert_eq!(input, output);
/// ```
pub fn exchange<T: PartialOrd>(arr: &mut [T]) {
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

/// Sort *arr* inplace using Bubble sort
pub fn bubble<T: PartialOrd>(arr: &mut [T]) {
    let n = arr.len();

    for i in 0..n {
        let mut swapped = false;
        for j in 0..n - i - 1 {
            if arr[j] > arr[j + 1] {
                arr.swap(j, j + 1);
                swapped = true;
            }
        }
        if !swapped {
            break;
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

/// Cosort two arrays
pub fn co_sort<O: PartialOrd, T>(arr: &mut [O], co_arr: &mut [T]) {
    for i in 1..arr.len() {
        let mut j = i;
        while j > 0 && arr[j] < arr[j - 1] {
            arr.swap(j, j - 1);
            co_arr.swap(j, j - 1);
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

/// Sort *arr* inplace using quicksort.
///
/// # Example
/// ```
/// use scialg::sort::quick;
///
/// let mut input = vec![3, 4, 1, 2];
/// let output = vec![1, 2, 3, 4];
///
/// quick(&mut input);
///
/// assert_eq!(input, output);
/// ```
pub fn quick<T: PartialOrd + Copy>(arr: &mut [T]) {
    let n = arr.len();
    let m: usize = 7;
    let mut il = 0;
    let mut ir = n - 1;
    let mut stack = Vec::new();

    loop {
        if ir - il < m {
            for i in 1..arr.len() {
                let mut j = i;
                while j > 0 && arr[j] < arr[j - 1] {
                    arr.swap(j, j - 1);
                    j -= 1;
                }
            }
            if stack.is_empty() {
                break;
            }
            ir = stack.pop().unwrap();
            il = stack.pop().unwrap();
            continue;
        }

        let k = (il + ir) / 2;
        arr.swap(k, il + 1);
        if arr[il] > arr[ir] {
            arr.swap(il, ir);
        }
        if arr[il + 1] > arr[ir] {
            arr.swap(il + 1, ir);
        }
        if arr[il] > arr[il + 1] {
            arr.swap(il, il + 1);
        }
        let mut i = il + 1;
        let mut j = ir;
        let pivot = arr[il + 1];
        loop {
            i += 1;
            j -= 1;
            while arr[i] < pivot {
                i += 1;
            }
            while arr[j] > pivot {
                j -= 1;
            }
            if j < i {
                break;
            }
            arr.swap(i, j);
        }
        arr[il + 1] = arr[j];
        arr[j] = pivot;

        if ir - i + 1 >= j - 1 {
            stack.push(i);
            stack.push(ir);
            ir = j - 1;
        } else {
            stack.push(il);
            stack.push(j - 1);
            il = i;
        }
    }
}
