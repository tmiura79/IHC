# adhe_a

This code calculates the measured Hubble constant based on the luminosity distance-redshift relation of a light source and used for our paper [1]. 
A perturbation with amplitude A and width L exists between the observer and the light source. 
The structure formation is calculated using the adhesion model.
Multiple calculations are performed with different values of A.
The quantities for each value of A are recorded in the file `output_aa.dat`.


## Requirements

- **C Compiler**: The code is written in C. The default compiler is `gcc`.

## Compilation and Execution

1. Compile the code using `make`:
    ```bash
    make
    ```
2. Run the program:
    ```bash
    ./adhe
    ```

## Licence

[MIT](https://github.com/tmiura79/IHC/blob/main/LICENSE)

## Reference
[1]　T. Miura and T. Tanaka, *JCAP* **05** (2024) 126, arXiv:2309.02288.

## Author

Taishi Miura

GitHub: [tmiura79](https://github.com/tmiura79)

