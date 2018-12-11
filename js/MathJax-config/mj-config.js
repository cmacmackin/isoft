window.MathJax = {
  AuthorInit: function() {
    MathJax.Hub.Register.StartupHook("TeX AMSmath Ready", function() {
      MathJax.InputJax.TeX.Definitions.Add({
        macros: {
          setcounter: "setCounter"
        }
      }, null, true);
      MathJax.InputJax.TeX.Parse.Augment({
        setCounter: function(name) {
          var num =  parseInt(this.GetArgument(name));
          MathJax.Extension["TeX/AMSmath"].number = num;
        }
      });
    });
  }
};
