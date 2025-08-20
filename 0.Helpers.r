Helper1 <- function(INPUT, UNIT_ID, TERM_A) {
  INPUT %>% 
    dplyr::rename(Unit_ID = !!rlang::sym(UNIT_ID), 
                  Term_A = !!rlang::sym(TERM_A)) %>% 
    dplyr::count(Unit_ID, Term_A) %>% 
    dplyr::group_by(Unit_ID) %>% 
    dplyr::mutate(Value = n / sum(n)) %>% 
    dplyr::ungroup() %>% 
    tidyr::complete(Unit_ID, Term_A, fill = list(Value = 0)) %>% 
    dplyr::transmute(!!rlang::sym(UNIT_ID) := Unit_ID, 
                     FeatureGroup = paste0('Prop', TERM_A), 
                     Term_A, 
                     Feature = Term_A,
                     Value) %>% 
    
    dplyr::group_by(Feature) %>% 
    dplyr::filter(!all(is.na(Value))) %>% 
    dplyr::ungroup()
}

Helper2 <- function(INPUT, UNIT_ID, TERM_A, TERM_B, PREPOSITION = 'in', 
                    REL2A, REL2B) {
  if (REL2A && REL2B) {
    stop("Both REL2A and REL2B cannot be TRUE.")
  }
  if (REL2A | REL2B) {
    RELATIVE_JOIN <- INPUT %>%
      Helper1(UNIT_ID, TERM_A = ifelse(REL2A, TERM_A, TERM_B)) %>%
      dplyr::transmute(!!rlang::sym(UNIT_ID),
                       RelativeTo = Feature,
                       Prop4Filter = Value)
  }
  
  INPUT %>% 
    dplyr::rename(Unit_ID = !!rlang::sym(UNIT_ID), 
                  Term_A = !!rlang::sym(TERM_A), 
                  Term_B = !!rlang::sym(TERM_B)) %>% 
    dplyr::select(Unit_ID, Term_A, Term_B) %>%
    dplyr::count(Unit_ID, Term_A, Term_B) %>% 
    {if (REL2A && !REL2B) dplyr::group_by(., Unit_ID, Term_A) 
      else if (REL2B && !REL2A) dplyr::group_by(., Unit_ID, Term_B) 
      else dplyr::group_by(., Unit_ID)} %>% 
    dplyr::mutate(Value = n / sum(n)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct() %>% 
    tidyr::complete(Unit_ID, Term_A, Term_B, fill = list(Value = 0)) %>% 
    dplyr::transmute(!!rlang::sym(UNIT_ID) := Unit_ID, 
                     FeatureGroup = paste0('Prop', TERM_A, 'in', TERM_B) %>% 
                       {if (REL2A) paste0(., 'rel2', TERM_A) 
                         else if (REL2B) paste0(., 'rel2', TERM_B)
                         else .}, 
                     Term_A, 
                     Term_B, 
                     Feature = paste0(Term_A, ';', Term_B) %>% 
                       {if (REL2A) paste0(., ';rel2_', Term_A) 
                         else if (REL2B) paste0(., ';rel2_', Term_B)
                         else .}, 
                     Value) %>% 
    {if (REL2A | REL2B) {
      dplyr::mutate(., RelativeTo = if (REL2A) Term_A else Term_B) %>% 
        dplyr::left_join(RELATIVE_JOIN, by = c(UNIT_ID, 'RelativeTo'))
    } else .} %>%
    dplyr::group_by(Feature) %>% 
    dplyr::filter(!all(is.na(Value))) %>% 
    dplyr::ungroup() %>% 
    dplyr::relocate(any_of(c(UNIT_ID, 'FeatureGroup', 'Term_A', 'Term_B', 'Attributive', 'Feature', 'Value', 'RelativeTo', 'Prop4Filter'))) %>% 
    dplyr::group_by(Feature) %>% 
    dplyr::mutate(ndis = dplyr::n_distinct(Value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(ndis > 1) %>% 
    dplyr::select(-ndis)
}
