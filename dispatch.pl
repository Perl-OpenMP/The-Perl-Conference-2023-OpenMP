
# caller
deux { foo => \&foo }
deux { bar => \&bar }
deux { baz => \&baz }
when {
  my $v = shift;
  if ($v < 0) {
    return q{foo};
  elsif ($v >= 0 and $v < 10) {
    return q{bar};
  else {
    return q{baz};
  } 
}
call {
  my $subref = shift;
  $subref->();
}, $value;
