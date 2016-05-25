<?php

namespace vladkolodka\linearAlgebra;

class MatrixException extends \Exception {
    private $errors = [
        'POLYMORPHIC', 'ARGUMENT_TYPE', 'ARGUMENT_BOUNDS', 'MATRIX_DIMENSION', 
        'ARRAY_LENGTH', 'RANK', 'MATRIX_SINGULAR', 'MATRIX_SQUARE', 'Matrix_SPD'
    ];

    public function __construct($code) {
        parent::__construct($this->errors[$code], $code, $this);
    }
}