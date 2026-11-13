// Left/right arrow keys go to the previous/next page — the targets of the
// previous/next cards at the foot of every page (overrides/partials/
// mc-page-footer.html; the home page hides its cards but keeps them in the
// markup for this). Clicking the card lets instant navigation take over.
//
// The keys only act when nothing that uses them has focus: the page itself,
// or a plain link or button. A focused field, code block (they scroll
// sideways), content tab, or the open search dialog (focus sits in its
// shadow host) keeps its arrows.
(() => {
  if (window.__mcArrowNavigation) return;
  window.__mcArrowNavigation = true;

  const keeps_arrows = "pre, .tabbed-set, [role=tablist], [role=listbox], [role=menu], [role=slider]";

  const arrowsAreFree = () => {
    const el = document.activeElement;
    if (!el || el === document.body || el === document.documentElement) return true;
    return el.matches("a[href], button") && !el.shadowRoot && !el.closest(keeps_arrows);
  };

  document.addEventListener("keydown", (event) => {
    if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
    if (event.defaultPrevented || event.altKey || event.ctrlKey || event.metaKey || event.shiftKey) return;
    if (!arrowsAreFree()) return;
    const direction = event.key === "ArrowLeft" ? "prev" : "next";
    const link = document.querySelector(`.mc-pagination__link--${direction}`);
    if (!link) return;
    event.preventDefault();
    link.click();
  });
})();
