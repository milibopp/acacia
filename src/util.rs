use nalgebra::{POrd, Orig};


/// Find infimum and supremum of a set of points
///
/// The function returns a tuple of the infimum and supremum. They are
/// determined component-wise according to the `POrd` implementation of the
/// point type.
///
/// ```
/// # extern crate nalgebra;
/// # extern crate acacia;
/// # fn main() {
/// use nalgebra::Pnt2;
/// use acacia::util::limits;
/// let points = vec![Pnt2::new(0.0f32, 1.0), Pnt2::new(2.0, -3.0)];
/// let (inf, sup) = limits(points.iter());
/// # }
/// ```
///
/// # Parameters
///
/// - `values` is a by-reference iterator over points
///
/// TODO: maybe this can be integrated into nalgebra?
pub fn limits<'a, I, P>(values: I) -> (P, P)
    where I: Iterator<&'a P>,
          P: POrd + Orig + Copy + 'a,
{
    let mut values = values;
    match values.next() {
        Some(&first) => values.fold(
            (first, first),
            |(inf, sup), new|
                (POrd::inf(&inf, new), POrd::sup(&sup, new))
        ),
        None => (Orig::orig(), Orig::orig()),
    }
}


#[cfg(test)]
mod test {
    use super::limits;
    use nalgebra::{Orig, Pnt3, POrd};

    #[test]
    fn limits_no_points() {
        let (inf, sup): (Pnt3<f64>, Pnt3<f64>) = limits(vec![].iter());
        assert_eq!(inf, Orig::orig());
        assert_eq!(sup, Orig::orig());
    }

    #[test]
    fn limits_one_point() {
        let p = Pnt3::new(1.0f64, 2.0, -3.0);
        let (inf, sup) = limits(vec![p].iter());
        assert_eq!(inf, p);
        assert_eq!(sup, p);
    }

    #[quickcheck]
    fn limits_inf_less_than_or_equal_sup(points: Vec<(f32, f32, f32)>) -> bool {
        let points: Vec<Pnt3<f32>> = points.iter().map(|&(x, y, z)| Pnt3::new(x, y, z)).collect();
        let (inf, sup) = limits(points.iter());
        POrd::partial_le(&inf, &sup)
    }

    #[quickcheck]
    fn limits_inf_less_than_or_equal_all(points: Vec<(f32, f32, f32)>) -> bool {
        let points: Vec<Pnt3<f32>> = points.iter().map(|&(x, y, z)| Pnt3::new(x, y, z)).collect();
        let (inf, _) = limits(points.iter());
        points.iter().all(|p| POrd::partial_le(&inf, p))
    }

    #[quickcheck]
    fn limits_sup_greater_than_or_equal_all(points: Vec<(f32, f32, f32)>) -> bool {
        let points: Vec<Pnt3<f32>> = points.iter().map(|&(x, y, z)| Pnt3::new(x, y, z)).collect();
        let (_, sup) = limits(points.iter());
        points.iter().all(|p| POrd::partial_ge(&sup, p))
    }
}
