use std::{env, fs, path::Path};

fn collect(root: &Path, directory: &Path, names: &mut Vec<String>) {
    for entry in fs::read_dir(directory).expect("read embedded viewer assets") {
        let path = entry.expect("viewer asset entry").path();
        if path.is_dir() {
            collect(root, &path, names);
        } else {
            names.push(
                path.strip_prefix(root)
                    .unwrap()
                    .to_string_lossy()
                    .replace('\\', "/"),
            );
        }
    }
}

fn main() {
    println!("cargo:rerun-if-changed=assets/viewer");
    let root = Path::new("assets/viewer");
    let mut names = Vec::new();
    collect(root, root, &mut names);
    names.sort();
    let mut source = String::from("static ASSETS: &[(&str, &[u8])] = &[\n");
    for name in names {
        source.push_str(&format!(
            "({name:?}, include_bytes!(concat!(env!(\"CARGO_MANIFEST_DIR\"), \"/assets/viewer/\", {name:?}))),\n"
        ));
    }
    source.push_str("];\n");
    fs::write(
        Path::new(&env::var_os("OUT_DIR").unwrap()).join("viewer-assets.rs"),
        source,
    )
    .expect("write embedded viewer index");
}
