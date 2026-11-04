// Solve a small LP and MIP through the highs crate, linked against an
// installed HiGHS found by pkg-config (highs-sys "discover" feature)
use highs::{HighsModelStatus, RowProblem, Sense};
use std::ffi::CStr;

fn solve(integer_z: bool) -> f64 {
    let mut pb = RowProblem::default();
    let x = pb.add_column(1., 0.0..=4.5);
    let y = pb.add_column(4., 1.0..=1.0);
    let z = if integer_z {
        pb.add_integer_column(9., 0.0..)
    } else {
        pb.add_column(9., 0.0..)
    };
    pb.add_row(..=5.5, &[(x, 1.), (y, 1.)]);
    pb.add_row(10.0.., &[(x, 1.), (z, 1.)]);
    pb.add_row(..=7.5, &[(y, -1.), (z, 1.)]);

    let mut model = pb.optimise(Sense::Minimise);
    model.make_quiet();
    let solved = model.solve();
    assert_eq!(solved.status(), HighsModelStatus::Optimal);
    solved.objective_value()
}

fn main() {
    let githash = unsafe { CStr::from_ptr(highs_sys::Highs_githash()) }
        .to_str()
        .unwrap();
    println!("Linked HiGHS git hash: {githash}");
    if let Ok(expected) = std::env::var("HIGHS_EXPECTED_GITHASH") {
        assert!(
            expected.starts_with(githash) || githash.starts_with(&expected),
            "linked HiGHS {githash} is not the checkout {expected}"
        );
    }

    let lp = solve(false);
    let mip = solve(true);
    println!("LP obj = {lp}, MIP obj = {mip}");
    assert!((lp - 58.).abs() < 1e-6);
    assert!((mip - 62.).abs() < 1e-6);
}
