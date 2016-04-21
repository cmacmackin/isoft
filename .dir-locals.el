((nil . ((eval . (set (make-local-variable 'my-project-path)
                      (file-name-directory
                       (let ((d (dir-locals-find-file ".")))
                         (if (stringp d) d (car d))))))
         (eval . (message "Project directory set to `%s'." my-project-path))
         (eval . (setq flycheck-gfortran-include-path
              (mapcar (lambda (p) (expand-file-name p my-project-path)) 
                          '("~/.local/include/" 
	                    "./mod"
		            "./factual/mod"
			    "./tests/mod"))))
         (eval . (message "Include directories set to `%s'." flycheck-gfortran-include-path))
	 (eval . (setq flycheck-gfortran-args 
		     (concat "-J" 
		         (expand-file-name '"./mod" my-project-path))))
         (eval . (message "Gfortran arguments set to `%s'." flycheck-gfortran-args))
)))
